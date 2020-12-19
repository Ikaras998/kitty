/* kitty: C++ truth table library
 * Copyright (C) 2017-2020  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file threshold_identification.hpp
  \brief Threshold logic function identification

  \author CS-472 2020 Fall students
*/

#pragma once

#include <vector>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"
#include "isop.hpp"

namespace kitty
{

/*! \brief Threshold logic function identification

  Given a truth table, this function determines whether it is a threshold logic function (TF)
  and finds a linear form if it is. A Boolean function is a TF if it can be expressed as

  f(x_1, ..., x_n) = \sum_{i=1}^n w_i x_i >= T

  where w_i are the weight values and T is the threshold value.
  The linear form of a TF is the vector [w_1, ..., w_n; T].

  \param tt The truth table
  \param plf Pointer to a vector that will hold a linear form of `tt` if it is a TF.
             The linear form has `tt.num_vars()` weight values and the threshold value
             in the end.
  \return `true` if `tt` is a TF; `false` if `tt` is a non-TF.
*/
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
  std::vector<int64_t> linear_form;
  std::vector<bool> unateness;

  auto ttCopy = tt;


  //Checks if the function is binate, positive/ negative unate for each variable
  int numVars = tt.num_vars();
  int size = numVars*2;
  for(uint8_t i = 0; i < numVars; i++){

      auto cof1 = cofactor1( ttCopy, i );
      auto cof0 = cofactor0(ttCopy, i);

      int posUn =0 ;
      int negUn = 0;

      auto it2 = cof0.begin();
      for(auto it1 = cof1.begin(); it1 != cof1.end(); it1++, it2++){
          uint64_t limit= 1;

          for(int k = 0; k < size; k++){
              limit *= 2;
          }
          uint64_t cof1Bits = *it1;
          uint64_t cof0Bits = *it2;


          for(uint64_t j = 1; j <= limit; j*=2 ){
              uint64_t maskedCof1 = cof1Bits & j;
              uint64_t maskedCof0 = cof0Bits & j;


              if(maskedCof0 == 0 && maskedCof1 != 0){
                  posUn = 1;
              }

              if(maskedCof1 == 0 && maskedCof0 != 0){
                  negUn = 1;
              }

              if(posUn == 1 and negUn == 1){

                  return false;
              }

          }

      }
      if(posUn == 1){
          unateness.push_back(true);
      }
      if(negUn == 1){
          unateness.push_back(false);
      }

  }

  //if unateness has a 0, flip the corresponding variable
  for(int i = 0; i < numVars; i ++){
      if(unateness[i] == 0){
          flip_inplace( ttCopy, i);
      }
  }


  auto on = isop(ttCopy);
  auto off = isop(~ttCopy);

  lprec *lp;

  lp = make_lp(0,0);
  if(lp == NULL){
      fprintf(stderr, "Unable to create LP model\n");
      return 0;
  }

  delete_lp(lp);

  /* if tt is TF: */
  /* push the weight and threshold values into `linear_form` */
  if ( plf )
  {
    *plf = linear_form;
  }
  return true;
}



} /* namespace kitty */
