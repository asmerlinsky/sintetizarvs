/*
ht.c		Joe Mietus		Dec 19 2001

-------------------------------------------------------------------------------
ht: Perform Hilbert transform
Copyright (C) 2002 Joe Mietus

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the author by e-mail (joe at physionet dot org) or postal mail
(MIT Room E25-505A, Cambridge, MA 02139 USA).  For updates to this software,
please visit PhysioNet (http://www.physionet.org/).
_______________________________________________________________________________
Perform Hilbert transform and calculate instantaneous magnitude and frequency.
*/


#define LMAXIMO 560000
#define LFILT 128	/* impulse response filter length; must be even!!! */


hilbert(double *x, double *transformed, int size)
{

    int i, npt, lfilt;
    //static double time[LMAXIMO+1], x[LMAXIMO+1];
    static double *xh, *ampl, *phase, *omega;
    double hilb[LFILT+1], pi, pi2, xt, xht;

    xh=dvector(1,LMAXIMO);
    ampl=dvector(1,LMAXIMO);
    phase=dvector(1,LMAXIMO);
    omega=dvector(1,LMAXIMO);

    pi = M_PI;
    pi2 = 2*pi;
    lfilt = LFILT;

        

    for (i=1; i<=lfilt; i++)
        hilb[i]=1/((i-lfilt/2)-0.5)/pi;

    for (i=1; i<=size; i++) {
//    for (i=1; i<=LMAXIMO && scanf("%lf %lf", &time[i], &x[i])==2; i++) {
        xh[i] = 0.0;
        ampl[i] = 0.0;
    }
    npt=i-1;

    if (npt < lfilt){ 
	printf("insufficient data\n");
	exit(0);
    }

    /* Hilbert transform */
      convol(x, xh, hilb, npt, lfilt);

    /* Ampl and phase */
    for (i=lfilt/2+1; i<=npt-lfilt/2; i++) {
        xt = x[i];
        xht = xh[i];
        ampl[i] = sqrt(xt*xt+xht*xht);
        phase[i] = atan2(xht ,xt);
        if (phase[i] < phase[i-1])
            omega[i] = phase[i]-phase[i-1]+pi2;
        else
            omega[i] = phase[i]-phase[i-1];
    }

    for (i=lfilt/2+2; i<=npt-lfilt/2; i++) transformed[i]= ampl[i]; 
	//fase: omega[i]/pi2);

}


convol(source, target, filt, npt, lfilt)
int npt, lfilt;
double *source, *target, *filt;
{
    int i, l;
    double yt;

    for (l=1; l<=npt-lfilt+1; l++) {
        yt = 0.0;
        for (i=1; i<=lfilt; i++) 
            yt = yt+source[l+i-1]*filt[lfilt+1-i];
        target[l] = yt;
    }

    /* shifting lfilt/1+1/2 points */
    for (i=1; i<=npt-lfilt; i++) {
        target[i] = 0.5*(target[i]+target[i+1]);
    }
    for (i=npt-lfilt; i>=1; i--)
        target[i+lfilt/2]=target[i];

    /* writing zeros */
    for (i=1; i<=lfilt/2; i++) {
        target[i] = 0.0;
        target[npt+1-i] = 0.0;
    }
}



