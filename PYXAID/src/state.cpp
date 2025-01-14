/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
 ***********************************************************/

#include "state.h"
#include "namd.h"
#include "aux.h"

int ext2int(int external, vector<int> &active_space)
{
  // Orbital indexing conversion
  // This definition is not general
  // External orbitals: 1, -1, 2, -2, ... (1-is the lowest, + = alpha, - = betha)
  // Internal orbitals: 0,  1, 2,  3, ...
  // General: External orbital = ext_orb
  //          Internal orbital = int_orb
  //          related as: active_space[k] = |ext_orb|;
  //          and then: ext_orb = sign(ext_orb) * k
  //          and finally not general scheme is applied
  //
  int internal;
  for (int k = 0; k < active_space.size(); k++)
  {
    if (active_space[k] == abs(external))
    {
      internal = k;
      break;
    }
  }
  int f = (external > 0) ? 0 : 1;
  internal = 2 * abs(internal) + f;
  return internal;
}

int delta(vector<int> &A, vector<int> &B, int &a, int &b)
{
  // This is more general approach - search for pair of indices in A and B
  // which are different, they are then assigned to a and b correspondingly
  // All other indexes in A and B should be equivalent (not necessarily in the
  // same order, this is just set operations)
  // Returns 1 if only 1 pair of indices is different between A and B,
  // return 0 otherwise
  // That is res = 1 means A and B are coupled, while res = 0 - not
  // Algorithm: first construct the list of overlapping orbitals
  // second - go through this list and determine the populations of
  // each orbital in this list in each of the configurations
  //

  // THIS VERSION REGARDS SIGN OF THE ORBITAL INDEX

  int res = 1;
  int sz = A.size();
  int nA = 0;             // number of elements in A which do not exist in B
  int nB = 0;             // number of elements in B which do not exist in A
  vector<int> _A, _B, _C; // modules of A and B and overlap C
  for (int i = 0; i < sz; i++)
  { // the size of A and B is assumed to be the same
    int mA = (A[i]);
    int mB = (B[i]);
    _A.push_back(mA);
    _B.push_back(mB);
    if (!is_in_vector(mA, _C))
    {
      _C.push_back(mA);
    }
    if (!is_in_vector(mB, _C))
    {
      _C.push_back(mB);
    }
  } // for i

  sz = _C.size();

  int nexc = 0; // Number of excitations between 2 states
  for (int i = 0; i < sz; i++)
  {
    int n_in_a, n_in_b;
    vector<int> tmpa, tmpb;
    n_in_a = num_in_vector(_C[i], _A, tmpa);
    n_in_b = num_in_vector(_C[i], _B, tmpb);
    int d = n_in_a - n_in_b;
    if (d > 0)
    {
      nexc += d;
    }
    if (d == 1)
    {
      a = _C[i];
    }
    if (d == -1)
    {
      b = _C[i];
    }
  } // for i

  if (nexc == 1)
  {
    res = 1;
    // The orbitals should already be known, because only one pair is different
  }
  else
  {
    res = 0;
  }

  return res;
}

int me_state::calculate_Exc(vector<int> &Exc_i, vector<int> &Exc_j, vector<double> &Exc_val, vector<int> &shift_i, vector<double> &shift_E)
{
  Exc = 0.0;
  int Nel = actual_state.size();
  int sz = Exc_val.size(); //

  // Calculate correction due to the presence of more than 1 electron on given orbital
  for (int i = 0; i < Nel - 1; i++)
  {
    for (int j = i + 1; j < Nel; j++)
    {
      for (int n = 0; n < sz; n++)
      {
        if (((abs(actual_state[i]) == Exc_i[n]) && (abs(actual_state[j]) == Exc_j[n])) ||
            ((abs(actual_state[i]) == Exc_j[n]) && (abs(actual_state[j]) == Exc_i[n])))
        {
          Exc += Exc_val[n];
        }
      } // for n
    }   // for j
  }     // for i

  // Shift given (1-electron orbitals) - "Scissor" operator
  sz = shift_E.size();
  for (int i = 0; i < Nel; i++)
  {
    for (int n = 0; n < sz; n++)
    {
      if (abs(actual_state[i]) == shift_i[n])
      {
        Exc += shift_E[n];
      }
    } // for n
  }   // for i

  return 1;
}

void me_state::show_state()
{
  cout << "Active space: ";
  show_vector(active_space);
  cout << endl;
  cout << "Actual state: ";
  show_vector(actual_state);
  cout << endl;
}

int list2state(json lst, vector<int> &active_space, me_state &ES)
{
  //        string   int list   double(optinal)
  // lst = [  name,  [1,2,3],     Eshift   ]
  int res = 1;
  int sz = lst.size();
  if (sz >= 2)
  {
    // Name - field 0
    ES.name = lst[0].get<std::string>();

    // Actual state - field 1
    vector<int> lst1 = lst[1].get<vector<int>>();
    int sz1 = lst1.size();
    vector<int> state;
    for (int i = 0; i < sz1; i++)
    {
      if (is_in_vector(abs(lst1[i]), active_space))
      {
        state.push_back(lst1[i]);
      }
      else
      {
        res = 0;
        break;
      }
    } // for i
    if (res)
    {
      ES.set_me_state(active_space, state);
    }

    // Eshift - field 2 (optional)
    if (sz >= 3)
    {
      ES.Eshift = lst[2].get<double>();
    }
  }
  else
  {
    cout << "Format Error(in list2state): the state is given as a list of at least 2 entries:";
    cout << "[<label>(string), <actual state>(list of integers), <Eshift>(float, optional) ]\n";
  }

  return res;
}

void input_states(json params, vector<me_state> &states)
{
  // States are defined as a list in a dictionary with a key "states"
  // For example:
  // params["states"] = []
  // params["states"].append(["GS",[12,-12,13,-13], 0.00])      # 0
  //
  // params["states"].append(["S1",[12,-12,-13,19],S1_corr])    # 1
  // params["states"].append(["S1a",[-12,13,-13,19],S1_corr])   # 2

  int is_active_space, is_ground_state, is_excl_ground_state;
  is_active_space = is_ground_state = is_excl_ground_state = 0; // Not yet defined
  me_state GS, ES;                                              // This is a ground state and excited state
  vector<int> active_space, ground_state;
  vector<int> Exc_i, Exc_j;         // indexes of orbitals for given correlation correction
  vector<double> Exc;               // correlation correction
  vector<int> shift_i;              // indexes of the (1-electron) orbitals to be shifted
  vector<double> shift_E;           // energy by which shift corresponding (1-electron) orbitals
  vector<int> nac_scl_i, nac_scl_j; // indexes of the macrostates
  vector<double> nac_scl;           // scaling constant for given pair of the macrostates
  // First - look only for active space
  if (!params["active_space"].is_null())
  {
    active_space = params["active_space"].get<vector<int>>();
    is_active_space = 1;
  }

  // Now read the microstates and create corresponding determinants
  if (!params["states"].is_null() && is_active_space)
  {
    json micro = params["states"];
    for (int j = 0; j < micro.size(); j++)
    {
      json tmp = micro[j];
      if (list2state(tmp, active_space, ES))
      {
        states.push_back(ES);
      }
    }
  }

  // Read other orbital/determinant parameters
  if (!params["shift"].is_null())
  {
    int sz = params["shift"].size();
    shift_i = vector<int>(sz, 0);
    shift_E = vector<double>(sz, 0.0);
    for (int i = 0; i < sz; i++)
    {
      shift_i[i] = params["shift"][i][0].get<int>();
      shift_E[i] = params["shift"][i][1].get<double>();
    }
  }
  if (!params["Exc"].is_null())
  {
    int sz = params["Exc"].size();
    Exc_i = vector<int>(sz, 0);
    Exc_j = vector<int>(sz, 0);
    Exc = vector<double>(sz, 0.0);
    for (int i = 0; i < sz; i++)
    {
      Exc_i[i] = params["Exc"][i][0].get<int>();
      Exc_j[i] = params["Exc"][i][1].get<int>();
      Exc[i] = params["Exc"][i][2].get<double>();
    }
  }
  if (!params["nac_scale"].is_null())
  {
    int sz = params["nac_scale"].size();
    nac_scl_i = vector<int>(sz, 0);
    nac_scl_j = vector<int>(sz, 0);
    nac_scl = vector<double>(sz, 0.0);
    for (int i = 0; i < sz; i++)
    {
      nac_scl_i[i] = params["nac_scale"][i][0].get<int>();
      nac_scl_j[i] = params["nac_scale"][i][1].get<int>();
      nac_scl[i] = params["nac_scale"][i][2].get<double>();
    }
  }
  // Now calculate the Exc corrections for all states, based on extracted parameters Exc_i, Exc_j and Exc
  int sz = states.size();
  for (int i = 0; i < sz; i++)
  {
    states[i].calculate_Exc(Exc_i, Exc_j, Exc, shift_i, shift_E);
  }

  // Debugging
  cout << "Number of basis multi-electron states is: " << sz << endl;
  for (int i = 0; i < sz; i++)
  {
    cout << "State " << i << " : ";
    states[i].show_state();
    cout << " Exc = " << states[i].Exc << endl;
  }

  // Now set NAC scalings
  int sz1 = nac_scl.size(); // total number of all pairs to be scales
  for (int i = 0; i < sz; i++)
  { // For each state i
    for (int k = 0; k < sz1; k++)
    { // check all pairs j
      if (nac_scl_i[k] == i)
      { // this will scale d[nac_scl_i[k]==i][nac_scl_j[k]]
        // Check if state nac_scl_j[k] is already in list
        if (is_in_vector(nac_scl_j[k], states[i].nac_scl_indx))
        {
        }
        else
        { // Not yet included
          states[i].nac_scl_indx.push_back(nac_scl_j[k]);
          states[i].nac_scl.push_back(nac_scl[k]);
        }
      }
      else if (nac_scl_j[k] == i)
      { // this will scale d[nac_scl_i[k]][nac_scl_j[k]==i]
        // Check if state nac_scl_j[k] is already in list
        if (is_in_vector(nac_scl_i[k], states[i].nac_scl_indx))
        {
        }
        else
        { // Not yet included
          states[i].nac_scl_indx.push_back(nac_scl_i[k]);
          states[i].nac_scl.push_back(nac_scl[k]);
        }
      }
      else
      {
      } // skip that pair
    }   // for k
  }     // for i

  // Now print all corrections:
  for (int i = 0; i < sz; i++)
  {
    cout << "Couplings of the macrostate " << i << " will be scaled for the following states:\n";
    for (int k = 0; k < states[i].nac_scl.size(); k++)
    {
      cout << "     " << states[i].nac_scl_indx[k] << "   " << states[i].nac_scl[k] << endl;
    }
  } // for i
}

void input_iconds(json params, int me_numstates, vector<vector<int>> &iconds)
{
  // initial condition are defined as a list in a dictionary with a key "iconds"

  iconds = params["iconds"].get<vector<vector<int>>>();

  for (int i = 0; i < iconds.size(); i++)
  {
    if (iconds[i][1] < 0)
    {
      cout << "Error: Minimal excitation state is 0 (0 - is a ground state)\n";
      exit(0);
    }
    if (iconds[i][1] > me_numstates)
    {
      cout << "Error: The initial excitation state must be in range [ 0 , " << me_numstates << ")\n";
      exit(0);
    }
  }
}
