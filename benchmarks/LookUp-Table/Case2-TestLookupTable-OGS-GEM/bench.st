GeoSys-ST: Source Terms ------------------------------------------------

#STOP

#SOURCE_TERM
 $PCS_TYPE
  MASS_TRANSPORT
 $PRIMARY_VARIABLE 
  3-C
 $GEO_TYPE
  POINT POINT0
 $DIS_TYPE
 CONSTANT   0.015e-4   ; 0.015 mol m-3 a-1 4.753215e-10 mol m-3 s-1
#SOURCE_TERM
 $PCS_TYPE
  MASS_TRANSPORT
 $PRIMARY_VARIABLE 
  12-O
 $GEO_TYPE
  POINT POINT0
 $DIS_TYPE
  CONSTANT  0.03e-4  ; aprox CO2 rate : 0.03 mol m-3 a-1 9.50643e-10 mol m-3 s-1
#STOP

#STOP

