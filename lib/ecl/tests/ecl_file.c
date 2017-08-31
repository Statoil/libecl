
/*
   Copyright (C) 2013  Statoil ASA, Norway.

   The file 'ecl_file.c' is part of ERT - Ensemble based Reservoir Tool.

   ERT is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ERT is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.

   See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
   for more details.
*/
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>

#include <ert/util/test_util.h>
#include <ert/util/util.h>
#include <ert/util/test_work_area.h>

#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_file_view.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_endian_flip.h>

void test_writable(size_t data_size) {
  test_work_area_type * work_area = test_work_area_alloc("ecl_file_writable");
  const char * data_file_name = "test_file";

  ecl_kw_type * kw = ecl_kw_alloc("TEST_KW", data_size, ECL_INT);
  for(int i = 0; i < data_size; ++i)
    ecl_kw_iset_int(kw, i, ((i*37)+11)%data_size);

  fortio_type * fortio = fortio_open_writer(data_file_name, false, true);
  ecl_kw_fwrite(kw, fortio); 
  fortio_fclose(fortio);

  for(int i = 0; i < 4; ++i) {
    ecl_file_type * ecl_file = ecl_file_open(data_file_name, ECL_FILE_WRITABLE);
    ecl_kw_type * loaded_kw = ecl_file_view_iget_kw(
                                  ecl_file_get_global_view(ecl_file),
                                  0);
    test_assert_true(ecl_kw_equal(kw, loaded_kw));

    ecl_file_save_kw(ecl_file, loaded_kw);
    ecl_file_close(ecl_file);
  }

  ecl_kw_free(kw);
  test_work_area_free( work_area );
}

void test_truncated() {
  test_work_area_type * work_area = test_work_area_alloc("ecl_file_truncated" );
  {
    ecl_grid_type * grid = ecl_grid_alloc_rectangular(20,20,20,1,1,1,NULL);
    ecl_grid_fwrite_EGRID2( grid , "TEST.EGRID", ECL_METRIC_UNITS );
    ecl_grid_free( grid );
  }
  {
    ecl_file_type * ecl_file = ecl_file_open("TEST.EGRID" , 0 );
    test_assert_true( ecl_file_is_instance( ecl_file ) );
    ecl_file_close( ecl_file );
  }

  {
    offset_type file_size = util_file_size( "TEST.EGRID");
    FILE * stream = util_fopen("TEST.EGRID" , "r+");
    util_ftruncate( stream , file_size / 2 );
    fclose( stream );
  }
  {
    ecl_file_type * ecl_file = ecl_file_open("TEST.EGRID" , 0 );
    test_assert_NULL( ecl_file );
  }
  test_work_area_free( work_area );
}

void test_return_copy() {
  // create some dummy kws
  // create a dummy file w the kws
  // open ecl_file in RETURN_COPY mode
  // extract kws from the file, using ecl_file_get_kw
  // close ecl_file thus destroying kws within ecl_file
  // assert correct values in kws
  test_work_area_type * work_area = test_work_area_alloc("ecl_file_RETURN_COPY_testing");
  {
    char * file_name = "data_file";

    //creating the data file
    size_t data_size = 10;
    ecl_kw_type * kw1 = ecl_kw_alloc("TEST1_KW", data_size, ECL_INT);
    for(int i = 0; i < data_size; ++i)
       ecl_kw_iset_int(kw1, i, 537 + i);
    fortio_type * fortio = fortio_open_writer(file_name, false, ECL_ENDIAN_FLIP);
    ecl_kw_fwrite(kw1, fortio); 
   
    data_size = 5;
    ecl_kw_type * kw2 = ecl_kw_alloc("TEST2_KW", data_size, ECL_FLOAT);
    for(int i = 0; i < data_size; ++i)
      ecl_kw_iset_float(kw2, i, 0.15 * i);
    ecl_kw_fwrite(kw2, fortio);
    fortio_fclose(fortio); 
    //finished creating data file

    ecl_file_type * ecl_file = ecl_file_open(file_name, ECL_FILE_RETURN_COPY);
    ecl_kw_type * kw1_copy = ecl_file_iget_kw( ecl_file , 0 );
    ecl_kw_type * kw2_copy = ecl_file_iget_kw( ecl_file , 1 );

    test_assert_NULL(kw1_copy);
    //printf(" ******** %s: kw_size: %d\n", ecl_kw_get_size( kw1_copy ));
    ecl_file_close( ecl_file ); 

    //test_assert_true (ecl_kw_equal(kw1, kw1_copy));
    //test_assert_true (ecl_kw_equal(kw2, kw2_copy));
    
  }
  test_work_area_free( work_area );

}

int main( int argc , char ** argv) {
  util_install_signals();
  //test_writable(10);
  //test_writable(1337);
  //test_truncated();
  test_return_copy();
  exit(0);
}
