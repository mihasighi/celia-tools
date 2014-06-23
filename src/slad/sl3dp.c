/**************************************************************************/
/*                                                                        */
/*  CINV Library / Shape Domain                                           */
/*                                                                        */
/*  Copyright (C) 2009-2011                                               */
/*    LIAFA (University of Paris Diderot and CNRS)                        */
/*                                                                        */
/*                                                                        */
/*  you can redistribute it and/or modify it under the terms of the GNU   */
/*  Lesser General Public License as published by the Free Software       */
/*  Foundation, version 3.                                                */
/*                                                                        */
/*  It is distributed in the hope that it will be useful,                 */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*  GNU Lesser General Public License for more details.                   */
/*                                                                        */
/*  See the GNU Lesser General Public License version 3.                  */
/*  for more details (enclosed in the file LICENSE).                      */
/*                                                                        */
/**************************************************************************/

#include <stdio.h>
#include "smtlib2sl3.h"
#include "shad.h"

/**
 * Print informations on usage.
 */
void
print_help (void)
{
  printf ("sl3dp: decision procedure for SL3, version 0.1\n");
  printf ("Usage: sl3dp <file>\n");
  printf ("\t<file>: input file in the SMTLIB2 format");
  printf ("See http://www.liafa.jussieu.fr/celia/ for more details.");
}

/**
 * Entry of the decision procedure.
 * Call: sl3dp file.smt2
 */
int
main (int argc, char** argv)
{
  // Step 0: Check the arguments
  if (argc <= 1)
    {
      print_help ();
      return 1;
    }
  // Parse the file and response to check-sat
  // pre: the file shall exists.
  FILE* f = fopen (argv[1], "r");
  if (!f)
    {
      printf ("File %s not found!\nquit.", argv[1]);
      return 1;
    }
  smtlib2_sl3_parser *sp = smtlib2_sl3_parser_new ();
  smtlib2_abstract_parser_parse ((smtlib2_abstract_parser *) sp, f);
  smtlib2_sl3_parser_delete (sp);
  fclose (f);
  // post: sh_pos an sh_neg are initialized with the formulas read.
  // response to check-sat
  return 0;
}
