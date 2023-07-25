/***********************************************************
 * Copyright (C) 2013 PYXAID group
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
 ***********************************************************/

#include <iostream>
#include "namd.h"
using namespace std;

void version()
{
  cout << "================================= Information ==================================\n";
  cout << "PYXAID: PYthon eXtension for Ab Inition Dynamics version 1.0\n";
  cout << "/***********************************************************\n";
  cout << " * Copyright (C) 2013 PYXAID group\n";
  cout << " * This program is free software distributed under the terms of the\n";
  cout << " * GNU General Public License as published by the\n";
  cout << " * Free Software Foundation; either version 3 of the\n";
  cout << " * License, or (at your option) any later version.\n";
  cout << " * http://www.gnu.org/copyleft/gpl.txt\n";
  cout << "***********************************************************/\n";
}

void developers()
{
  cout << "===== Name ===================== Developers ======== Contact info ==============\n";
  cout << "Alexey V. Akimov                             e-mail: alexvakimov@gmail.com      \n";
  cout << "Oleg V. Prezhdo                              e-mail: oprezhdo@chem.rochester.edu\n";
  cout << "Hsiaoyi Tsai                                 e-mail: tsaihsiaoyi0228@gmail.com  \n";
  cout << "...\n";
}

void documentation()
{
  cout << "================================ Documentation =================================\n";
  cout << "Coming soon...\n";
}

int main()
{
  version();
  developers();
  documentation();
  namd();
  return 0;
}