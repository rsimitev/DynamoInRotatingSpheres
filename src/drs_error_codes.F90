! Copyright  L. Silva (lacsilva@gmail.com), 2015
module drs_error_codes
   implicit none
   integer, parameter:: WARN_TIME_STEP_TOO_BIG          = -200
   integer, parameter:: WARN_RES_TOO_LOW_PH             = -120
   integer, parameter:: WARN_RES_TOO_LOW_TH             = -110
   integer, parameter:: WARN_RES_TOO_LOW_R              = -100
   integer, parameter:: ERR_NR1_NOT_POW2                =   10
   integer, parameter:: ERR_NT_TOO_SMALL                =   20
   integer, parameter:: ERR_INVALID_M_SYMMETRY          =   30
   integer, parameter:: ERR_INCOMPATIBLE_SYMMETRY       =   40
   integer, parameter:: ERR_ILLEGAL_NPROCS              =  100
   integer, parameter:: ERR_MPI_NOT_INITIALISED         =  101
   integer, parameter:: ERR_TOO_MANY_PROCS_FOR_M        =  110
   integer, parameter:: ERR_TOO_MANY_PROCS_FOR_TH       =  120
   integer, parameter:: ERR_UNKNOWN_CONFIG_MAGIC_NUMBER =  200
   integer, parameter:: ERR_UNKNOWN_TEMP_BC             =  201
   integer, parameter:: ERR_UNKNOWN_COMP_BC             =  202
   integer, parameter:: ERR_UNKNOWN_FLOW_BC             =  203
   integer, parameter:: ERR_UNKNOWN_MAG_BC              =  204
   integer, parameter:: ERR_UNKNOWN_CALC_TYPE           =  910
   integer, parameter:: ERR_CONFIG_UNOPENABLE           =  300
   integer, parameter:: ERR_UNKNOWN_VARIABLE            =  301
   integer, parameter:: ERR_NO_COMP_SUPPORT             =  302
   integer, parameter:: ERR_CREATING_LOCK               =  995
   integer, parameter:: ERR_DELETING_LOCK               =  996
   integer, parameter:: ERR_UNCONFIGURED_LOCK           =  997
   integer, parameter:: ERR_NO_PAR_FILE                 = 2001
   integer, parameter:: ERR_NO_FLOW_POL_FILE            = 2002
   integer, parameter:: ERR_NO_FLOW_TOR_FILE            = 2004
   integer, parameter:: ERR_NO_MAG_POL_FILE             = 2016
   integer, parameter:: ERR_NO_MAG_TOR_FILE             = 2032
   integer, parameter:: ERR_NO_TEMPERATURE_FILE         = 2008
   integer, parameter:: ERR_NO_COMPOSITION_FILE         = 2064
   integer, parameter:: ERR_RES_TOO_LOW_FOR_HARTMAN     = 5000
   integer, parameter:: ERR_KENERGY_IS_NAN              = 6000
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
