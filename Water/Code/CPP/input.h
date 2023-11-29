/*
*     Nested L-shaped method for DRSO and Risk-Averse Optimization
*
*     VERSION 0.1
*
*     Authors:   Hamed Rahimian
*                The Ohio State University
*
*        Guzin Bayraksan and Tito Homem-de-Mello
*
*		September 21, 2017
*
*     (C)opyright 2017 - H. Rahimian, G. Bayraksan, and T. Homem-de-Mello
*
*/

#ifndef INPUT_H
#define INPUT_H


#include "types.h"

void ReadData(int& numChild, const char progname[]);
void ReadCoeff(const char ambiguitytype[], int& numInstance, std::vector < std::pair<double, double> >& config);

char* setFilename(void);

void setTypes(char **argv);

VariantType getVariantType(void);

#endif
