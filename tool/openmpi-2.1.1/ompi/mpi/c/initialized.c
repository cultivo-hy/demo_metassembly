/*
 * Copyright (c) 2004-2007 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2005 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2015      Research Organization for Information Science
 *                         and Technology (RIST). All rights reserved.
 * Copyright (c) 2015      Cisco Systems, Inc.  All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "ompi_config.h"

#include "ompi/mpi/c/bindings.h"
#include "ompi/runtime/params.h"
#include "ompi/communicator/communicator.h"
#include "ompi/errhandler/errhandler.h"

#if OMPI_BUILD_MPI_PROFILING
#if OPAL_HAVE_WEAK_SYMBOLS
#pragma weak MPI_Initialized = PMPI_Initialized
#endif
#define MPI_Initialized PMPI_Initialized
#endif

static const char FUNC_NAME[] = "MPI_Initialized";


int MPI_Initialized(int *flag)
{
    MPI_Comm null = NULL;

    /* We must obtain the lock to guarnatee consistent values of
       ompi_mpi_initialized and ompi_mpi_finalized.  Note, too, that
       this lock is held for the bulk of the duration of
       ompi_mpi_init() and ompi_mpi_finalize(), so when we get the
       lock, we are guaranteed that some other thread is not part way
       through initialization or finalization. */
    opal_mutex_lock(&ompi_mpi_bootstrap_mutex);

    if (MPI_PARAM_CHECK) {
        if (NULL == flag) {

            /* If we have an error, the action that we take depends on
               whether we're currently (after MPI_Init and before
               MPI_Finalize) or not */

            if (ompi_mpi_initialized && !ompi_mpi_finalized) {
                opal_mutex_unlock(&ompi_mpi_bootstrap_mutex);
                return OMPI_ERRHANDLER_INVOKE(MPI_COMM_WORLD, MPI_ERR_ARG,
                                              FUNC_NAME);
            } else {
                opal_mutex_unlock(&ompi_mpi_bootstrap_mutex);
                return OMPI_ERRHANDLER_INVOKE(null, MPI_ERR_ARG,
                                              FUNC_NAME);
            }
        }
    }

    *flag = ompi_mpi_initialized;
    opal_mutex_unlock(&ompi_mpi_bootstrap_mutex);

    return MPI_SUCCESS;
}
