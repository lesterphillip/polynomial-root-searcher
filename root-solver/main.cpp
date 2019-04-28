#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <complex>

using namespace std;

//FUNCTION DECLARATIONS:
complex <long double> poly(long double polynomial[], int degrees,
                           complex <long double> x);
complex <long double> *DurandKerner(long double polynomial[],
                                    int deg);

int main ()
{
    //CHECK FOR AVAILABILITY OF FILE AND OPEN IT
    string inFileName;
    cout << "Filename: ";
    cin >> inFileName;
    fstream inFileStr(inFileName.c_str(), ios::in);
    
    //CHECK IF FILE IS VALID
    if(inFileStr.fail())
    {
        cerr << "Unable to open " << inFileName << endl;
        system("pause");
        exit(0);
    }
    
    //CHECK IF FILE IS EMPTY
    if (inFileStr.peek() == EOF)
    {
        cout << "The file is empty. Exiting..." << endl;
        system("pause");
        return (0);
    }
    
    //CHECK IF FILE HAS INTEGERS IN IT
    int degree = 0;
    long double* coef = NULL;
    if (isdigit(inFileStr.peek()))
    {
        inFileStr >> degree;
        coef = new long double[degree];
        cout << "Degree: " << degree << endl;
        
        for (int i=0; i<=degree; i++)
        {
            inFileStr >> coef[i];
        }
    }
    
    //CHECK IF FILE HAS NON-NUMERIC CHARACTERS IN IT
    else
    {
        cerr << "You inputted an invalid file. Enter a .txt file "
        << "containing numbers only." << endl;
        system("pause");
        exit(0);
    }
    
    //ARRANGE COEFFICIENTS AND DIVIDE ALL BY LEADING COEFFICIENT
    long double* coefarrange = NULL;
    coefarrange = new long double[degree];
    cout << "Polynomial:" << endl;
    cout.precision(5);
    for (int i=degree, j=0; i>=0; i--, j++)
    {
        coefarrange[j] = coef[i]/coef[degree];
        cout << fixed << coefarrange[j] << "x^" << i << endl;
    }
    
    //GET THE ROOTS USING DURAND-KERNER METHOD
    complex <long double> *roots = DurandKerner(coefarrange, degree);
    
    //displays the roots calculated
    for(int print=0; print<degree; print++)
    {
        cout << scientific << "Root " << print+1 << " is " << roots[print]
        << endl;
    }
    
    //display-check if the roots are actual roots
    cout << "Checking:" << endl;
    for(int check=0; check<degree; check++)
    {
        cout << fixed << "f(" << roots[check] << ") = "
        << poly(coefarrange, degree, roots[check]) <<  endl;
    }
    
    system("pause");
    return 0;
}

complex <long double> *DurandKerner(long double polynomial[], int deg)
{
    complex<long double> constant(0.4,0.9), reset(1,0);
    //constant 0.4+0.9i is used for the formula of the code,
    //reset is used to reset the arrays below to 1, essential to
    //ensure y*=(a-b) works properly
    
    complex<long double> *rootApprox = NULL;
    complex<long double> firstDenom[deg], secondDenom[deg];
    rootApprox = new complex <long double> [deg];
    fill_n(firstDenom, deg, 1); //replaces all the values of the array with 1
    fill_n(secondDenom, deg, 1); //replaces all the values of the array with 1
    
    for (int i=0; i<deg; i++)
    {
        rootApprox[i] = pow(constant, i);
    }
    
    //iterates for a certain number of times until the root...
    //converges to the real number
    
    for (int iterations=0; iterations<=10000; iterations++)
    {
        //for loop goes through all the roots available
        //highest degree states the number of roots needed
        for(int nthRoot = 0; nthRoot<deg; nthRoot++)
        {
            //for the initializer
            //computes for the part of pn-1 - qn-1 and so on..
            for(int firstHalf=1; firstHalf<=deg-1; firstHalf++)
            {
                //ensures that the array starts at 1+0i for the first computation
                if(firstHalf==1)
                {
                    firstDenom[nthRoot]=reset;
                }
                
                if(nthRoot<firstHalf)
                {
                    firstDenom[nthRoot] *= (rootApprox[nthRoot] -
                                            rootApprox [firstHalf]);
                }
            }
            
            //initializes the polynomial function
            //f(qn-1), f(pn-1) ... etc
            complex <long double> y = poly(polynomial, deg, rootApprox[nthRoot]);
            
            //computes for  (qn-1 - pn) , (rn-1 - pn), ... etc
            //starts at nthRoot > 0 because formula says so
            if(nthRoot>0)
            {
                for(int secondHalf=0; secondHalf<=deg-1; secondHalf++)
                {
                    //ensures that the array starts at 1+0i for the first computation
                    if(secondHalf==0)
                    {
                        secondDenom[nthRoot]=reset;
                    }
                    
                    if(nthRoot>secondHalf)
                    {
                        secondDenom[nthRoot] *= (rootApprox[nthRoot] -
                                                 rootApprox[secondHalf]);
                    }
                }
            }
            //calculates for the formula of durand-kerner itself
            rootApprox[nthRoot] = (rootApprox [nthRoot] -
                                   (y/(firstDenom[nthRoot]*secondDenom[nthRoot])));
        }
    }
    return rootApprox;
}

complex <long double> poly(long double polynomial[], int degrees,
                           complex <long double> x)
{
    complex <long double> y = 0;
    for (int i = 0; i <= degrees; i++)
    {
        y += polynomial[i]*pow(x,(degrees-i));
    }
    return y;
}

