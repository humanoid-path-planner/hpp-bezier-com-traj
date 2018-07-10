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

  int getType(const VectorXd & mib, const VectorXd & mab, const int i)
  {
      const double& mibV = mib(i);
      const double& mabV = mab(i);
      int type = GLP_FR;
      if(mibV > UNBOUNDED_DOWN && mabV < UNBOUNDED_UP)
          type = GLP_DB;
      else if(mibV > UNBOUNDED_DOWN)
          type = GLP_LO;
      else if(mabV < UNBOUNDED_UP)
          type = GLP_UP;
      return type;
  }

  int solveglpk(const VectorXd & g0,
                    const MatrixXd & CE,
                    const VectorXd & ce0,
                    const MatrixXd & CI,
                    const VectorXd & ci0,
                    solvers::Cref_vectorX minBounds,
                    solvers::Cref_vectorX maxBounds,
                    VectorXd& x,
                    double& cost)
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
      const int numEqConstraints = (int)(CE.rows());
      const int numIneqConstraints =(int)(CI.rows());
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

      //COLUMNS
      glp_add_cols(lp, xsize);
      VectorXd miB =  minBounds.size() > 0 ? minBounds : VectorXd::Ones(xsize) * UNBOUNDED_DOWN;
      VectorXd maB =  maxBounds.size() > 0 ? maxBounds : VectorXd::Ones(xsize) * UNBOUNDED_UP;
      for (int i=0; i < xsize; ++i, ++idcol )
      {
          glp_set_col_bnds(lp, idcol, getType(miB, maB, i), miB(i), maB(i));
          glp_set_obj_coef(lp, idcol, g0(i));
      }
      glp_load_matrix(lp,idConsMat-1,ia,ja,ar);

      int res = glp_simplex(lp, &opts);
      if(res == 0)
      {
          cost = glp_get_obj_val(lp);                    //obtains a computed value of the objective function
          idrow = 1;
          for (int i =0; i < xsize; ++i, ++idrow)
              x(i) = glp_get_col_prim(lp, idrow);
      }
      glp_delete_prob(lp);                        //calls the routine glp_delete_prob, which frees all the memory
      glp_free_env();
      return res;
  }

  } /* namespace solvers */
