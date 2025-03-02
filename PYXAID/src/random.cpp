/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
 ***********************************************************/

#include "random.h"

//============================================================
//              Uniform distribution

double uniform(double a, double b)
{

  double ksi = rand() / ((double)RAND_MAX);
  return (a + (b - a) * ksi);
}
double p_uniform(double a, double b)
{
  return (1.0 / (b - a));
}

//============================================================
//              Exponential distribution

double exponential(double lambda)
{
  double ksi = uniform(0.0, 1.0);
  ksi = -(1.0 / lambda) * log(1.0 - ksi);
  return ksi;
}

double p_exp(double x, double lambda)
{
  return lambda * exp(-lambda * x);
}

//============================================================
//               Normal distribution

double normal()
{
  // Algorithm FL, Forsythe, Ahrens-Dieter
  // This algorithm generates Normally-distributed random variable as given by
  // function p_normal()
  double *a;
  a = new double[33];
  double *t;
  t = new double[33];
  double *h;
  h = new double[33];
  double *d;
  d = new double[48];

  a[1] = 0.00000000000000;
  t[1] = 0.00076738283767;
  h[1] = 0.03920617164634;
  a[2] = 0.03917608550309;
  t[2] = 0.00230687039764;
  h[2] = 0.03932704963665;
  a[3] = 0.07841241273311;
  t[3] = 0.00386061844387;
  h[3] = 0.03950999486086;
  a[4] = 0.11776987457909;
  t[4] = 0.00543845406707;
  h[4] = 0.03975702679515;
  a[5] = 0.15731068461017;
  t[5] = 0.00705069876857;
  h[5] = 0.04007092772490;
  a[6] = 0.19709908429430;
  t[6] = 0.00870839582019;
  h[6] = 0.04045532602655;
  a[7] = 0.23720210932878;
  t[7] = 0.01042356984914;
  h[7] = 0.04091480886081;
  a[8] = 0.27769043982157;
  t[8] = 0.01220953194966;
  h[8] = 0.04145507115859;
  a[9] = 0.31863936396437;
  t[9] = 0.01408124734637;
  h[9] = 0.04208311051344;
  a[10] = 0.36012989178957;
  t[10] = 0.01605578804548;
  h[10] = 0.04280748137995;
  a[11] = 0.40225006532172;
  t[11] = 0.01815290075142;
  h[11] = 0.04363862733472;
  a[12] = 0.44509652498551;
  t[12] = 0.02039573175398;
  h[12] = 0.04458931789605;
  a[13] = 0.48877641111466;
  t[13] = 0.02281176732513;
  h[13] = 0.04567522779560;
  a[14] = 0.53340970624127;
  t[14] = 0.02543407332319;
  h[14] = 0.04691571371696;
  a[15] = 0.57913216225555;
  t[15] = 0.02830295595118;
  h[15] = 0.04833486978119;
  a[16] = 0.62609901234641;
  t[16] = 0.03146822492920;
  h[16] = 0.04996298427702;
  a[17] = 0.67448975019607;
  t[17] = 0.03499233438388;
  h[17] = 0.05183858644724;
  a[18] = 0.72451438349236;
  t[18] = 0.03895482964836;
  h[18] = 0.05401138183398;
  a[19] = 0.77642176114792;
  t[19] = 0.04345878381672;
  h[19] = 0.05654656186515;
  a[20] = 0.83051087820539;
  t[20] = 0.04864034918076;
  h[20] = 0.05953130423884;
  a[21] = 0.88714655901887;
  t[21] = 0.05468333844273;
  h[21] = 0.06308488965373;
  a[22] = 0.94678175630104;
  t[22] = 0.06184222395816;
  h[22] = 0.06737503494905;
  a[23] = 1.00999016924958;
  t[23] = 0.07047982761667;
  h[23] = 0.07264543556657;
  a[24] = 1.07751556704027;
  t[24] = 0.08113194985866;
  h[24] = 0.07926471414968;
  a[25] = 1.15034938037600;
  t[25] = 0.09462443534514;
  h[25] = 0.08781922325338;
  a[26] = 1.22985875921658;
  t[26] = 0.11230007889456;
  h[26] = 0.09930398323927;
  a[27] = 1.31801089730353;
  t[27] = 0.13649799954975;
  h[27] = 0.11555994154118;
  a[28] = 1.41779713799625;
  t[28] = 0.17168856004707;
  h[28] = 0.14043438342816;
  a[29] = 1.53412054435253;
  t[29] = 0.22762405488269;
  h[29] = 0.18361418337460;
  a[30] = 1.67593972277344;
  t[30] = 0.33049802776911;
  h[30] = 0.27900163464163;
  a[31] = 1.86273186742164;
  t[31] = 0.58470309390507;
  h[31] = 0.70104742502766;
  a[32] = 2.15387469406144;
  t[32] = 0.00000000000000;
  h[32] = 0.00000000000000; // --   ---

  d[1] = 0.67448975019607;
  d[17] = 0.15040938382813;
  d[33] = 0.10597677198479;
  d[2] = 0.47585963017993;
  d[18] = 0.14590257684509;
  d[34] = 0.10433484129317;
  d[3] = 0.38377116397654;
  d[19] = 0.14177003276856;
  d[35] = 0.10276601206127;
  d[4] = 0.32861132306910;
  d[20] = 0.13796317369537;
  d[36] = 0.10126505151402;
  d[5] = 0.29114282663980;
  d[21] = 0.13444176150074;
  d[37] = 0.09982723448906;
  d[6] = 0.26368432217502;
  d[22] = 0.13117215026483;
  d[38] = 0.09844828202068;
  d[7] = 0.24250845238097;
  d[23] = 0.12812596512583;
  d[39] = 0.09712430874765;
  d[8] = 0.22556744380930;
  d[24] = 0.12527909006226;
  d[40] = 0.09585177768776;
  d[9] = 0.21163416577204;
  d[25] = 0.12261088288608;
  d[41] = 0.09462746119186;
  d[10] = 0.19992426749317;
  d[26] = 0.12010355965651;
  d[42] = 0.09344840710526;
  d[11] = 0.18991075842246;
  d[27] = 0.11774170701949;
  d[43] = 0.09231190933664;
  d[12] = 0.18122518100691;
  d[28] = 0.11551189226063;
  d[44] = 0.09121548217294;
  d[13] = 0.17360140038056;
  d[29] = 0.11340234879117;
  d[45] = 0.09015683778986;
  d[14] = 0.16684190866667;
  d[30] = 0.11140272044119;
  d[46] = 0.08913386650005;
  d[15] = 0.16079672918053;
  d[31] = 0.10950385201710;
  d[47] = 0.08814461935364;
  d[16] = 0.15534971747692;
  d[32] = 0.10769761656476;

  // cout<<"============== In normal ================\n";
  double res;
  double u = uniform(-1.0, 1.0);
  double y, A, up, w, T;
  int s;
  if (u < 0.0)
  {
    s = -1;
  }
  else if (u == 0.0)
  {
    s = 0;
  }
  else
  {
    s = 1;
  }
  u = fabs(u);
  u = 32 * u;
  int i = floor(u);

  if (i == 0)
  { // Tail
    i = 6;
    A = a[31]; // here I used 31 instead of 32 as in paper - this fixes the center-tail gap
    do
    {
      u = 2 * u;
      A = A + d[i];
      i = i + 1;
    } while (u <= 1.0);

    u = u - 1.0;
  id3:
    w = u * d[i];
    T = (0.5 * w + A) * w;
    while (1)
    {
      up = uniform(0.0, 1.0);
      u = uniform(0.0, 1.0);

      if (up >= T)
      {
        break;
      }
      else if (up <= u)
      {
        u = uniform(0.0, 1.0);
        goto id3;
      }
      else
      {
        T = u;
      }
    };

    y = A + w;
    res = s * y;
  }

  else
  { // Center
    up = u - i;
    A = a[i];

    while (1)
    {
    id2:
      if (up >= t[i])
      { // <--- B
        w = (up - t[i]) * h[i];
        break;
      }
      else
      {
        u = uniform(0.0, 1.0);
        w = u * (a[i + 1] - A);
        T = (0.5 * w + A) * w;
      id1:
        if (up >= T)
        {
          break;
        } // <--- A
        else
        {
          u = uniform(0.0, 1.0);

          if (up >= u)
          {
            T = u;
            up = uniform(0.0, 1.0);
            goto id1;
          } // goto A
          else
          {
            up = uniform(0.0, 1.0);
            goto id2;
          } // goto B
        }   // else
      }     // else
    }       // while

    y = A + w;
    res = s * y;

  } // center

  delete[] a;
  delete[] t;
  delete[] h;
  delete[] d;
  return res;
}
