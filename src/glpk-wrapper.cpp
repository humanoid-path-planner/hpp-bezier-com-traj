//
// Copyright (c) 2017 CNRS
//
// This file is part of tsid
// tsid is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
// tsid is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// tsid If not, see
// <http://www.gnu.org/licenses/>.
//


#include "solver/glpk-wrapper.hpp"
#include <glpk.h>
#include <iostream>


  namespace solvers
  {

  glpk_status solve(const VectorXd & g0,
                    const MatrixXd & CE,
                    const VectorXd & ce0,
                    const MatrixXd & CI,
                    const VectorXd & ci0,
                    VectorXd& x)
  {
      glp_smcp opts; glp_init_smcp(&opts); opts.msg_lev = GLP_MSG_OFF;
      glp_prob *lp;
      int ia[1 + 20000];                    //Row indices of each element
      int ja[1 + 20000];                    //column indices of each element
      double ar[1 + 20000];                     //numerical values of corresponding elements
      lp = glp_create_prob();                     //creates a problem object
      glp_set_prob_name(lp, "sample");            //assigns a symbolic name to the problem object
      glp_set_obj_dir(lp, GLP_MIN);               //calls the routine glp_set_obj_dir to set the
                                                  //omptimization direction flag,
                                                  //where GLP_MAX means maximization

      //ROWS
      // TODO SPECIFIC
      const int numEqConstraints = (int)(CE.rows());
      // for inequality, x.size()-3 last constraints are not relevant as they are the positivity constraints
      const int numIneqConstraints =(int)(CI.rows() - (x.size() -3));
      const int num_constraints_total (numEqConstraints + numIneqConstraints);
      glp_add_rows(lp, num_constraints_total);
      int idrow = 1;//adds three rows to the problem object
      int idcol = 1;
      int idConsMat = 1;
      int xsize = (int)(x.size());
      for (int i =0; i < numIneqConstraints; ++i, ++idrow )
      {
          glp_set_row_bnds(lp, idrow, GLP_UP, 0.0, ci0(i));
          for (int j =0; j < xsize; ++j, ++idcol )
          {
              if(CI(i,j) != 0.)
              {
                  ia[idConsMat] = idrow, ja[idConsMat] = idcol, ar[idConsMat] = CI(i,j); /* a[1,1] = 1 */
                  ++ idConsMat;
              }
          }
          idcol = 1;
      }
      for (int i =0; i < numEqConstraints; ++i, ++idrow )
      {
          glp_set_row_bnds(lp, idrow, GLP_FX, ce0(i), ce0(i));
          for (int j =0; j < xsize; ++j, ++idcol )
          {
              if(CE(i,j) != 0.)
              {
                  ia[idConsMat] = idrow, ja[idConsMat] = idcol, ar[idConsMat] = CE(i,j); /* a[1,1] = 1 */
                  ++ idConsMat;
              }
          }
          idcol = 1;
      }
      //assert(num_constraints_total == idConsMat -1);



      //COLUMNS
      glp_add_cols(lp, xsize); //adds three columns to the problem object
      // TODO SPECIFIC
      for (int i =0; i < 3; ++i, ++idcol )
      {
          glp_set_col_bnds(lp, idcol, GLP_FR, 0.0, 0.0);
          glp_set_obj_coef(lp, idcol, g0(i));
      }
      for (int i =3; i < xsize; ++i, ++idcol )
      {
          glp_set_col_bnds(lp, idcol, GLP_LO, 0.0, 0.0);
          glp_set_obj_coef(lp, idcol, g0(i));
      }
      /*for (int i =3; i < x.size(); ++i, ++idcol )
          glp_set_col_bnds(lp, idrow, GLP_UP, 0.0, ci0(i));*/

      // Constraint matrix
      /*Eigen::MatrixXd constraints = Eigen::MatrixXd::Zero(1+num_constraints_total,1+x.size());
      constraints.block(1,1,numIneqConstraints,idcol) = CI.block(0,0,numIneqConstraints,idcol);
      constraints.block(1+numIneqConstraints,1,numEqConstraints,idcol) = CE;*/

      glp_load_matrix(lp,idConsMat-1,ia,ja,ar);

      //ia[1] = 1, ja[1] = 1, ar[1] = 1.0; /* a[1,1] = 1 */

      int res = glp_simplex(lp, &opts);                      //calls the routine glp_simplex
      glpk_status ret = glpk_INFEASIBLE;
      if(res == 0)
      {
          ret = glpk_OPTIMAL;
                                                  //to solve LP problem
          //z = glp_get_obj_val(lp);                    //obtains a computed value of the objective function

          idrow = 1;
          for (int i =0; i < xsize; ++i, ++idrow)
              x(i) = glp_get_col_prim(lp, idrow);


          //printf("\nz = %g; x1 = %g; x2 = %g;\n", z, x1, x2); //writes out the optimal solution
      }
      glp_delete_prob(lp);                        //calls the routine glp_delete_prob, which frees all the memory
      glp_free_env();
      return ret;
  }

  } /* namespace solvers */
