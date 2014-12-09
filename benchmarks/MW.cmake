ADD_BENCHMARK ("MW" H_us/Darcy/unconf_WO_rch/unconf "OGS_FEM" 1
	H_us/Darcy/unconf_WO_rch/unconf_2.vtu
	H_us/Darcy/unconf_WO_rch/unconf_ply_bottom_t1.tec)

ADD_BENCHMARK ("MW" H_us/Darcy/unconf_W_rch/unconf "OGS_FEM" 1
	H_us/Darcy/unconf_W_rch/unconf_2.vtu
	H_us/Darcy/unconf_W_rch/unconf_ply_bottom_t1.tec)

ADD_BENCHMARK ("MW_LONG" DENSITY-DEPENDENT_FLOW/goswami-clement/constrBC_PressAsHead_tri/HM "OGS_FEM" 1
	DENSITY-DEPENDENT_FLOW/goswami-clement/constrBC_PressAsHead_tri/HM_46.vtu
	DENSITY-DEPENDENT_FLOW/goswami-clement/constrBC_PressAsHead_tri/HM_129.vtu
	DENSITY-DEPENDENT_FLOW/goswami-clement/constrBC_PressAsHead_tri/HM_168.vtu)

ADD_BENCHMARK ("MW_LONG" DENSITY-DEPENDENT_FLOW/goswami-clement/wholeBC_PressAsHead_tri/HM "OGS_FEM" 1
	DENSITY-DEPENDENT_FLOW/goswami-clement/wholeBC_PressAsHead_tri/HM_42.vtu
	DENSITY-DEPENDENT_FLOW/goswami-clement/wholeBC_PressAsHead_tri/HM_119.vtu
	DENSITY-DEPENDENT_FLOW/goswami-clement/wholeBC_PressAsHead_tri/HM_159.vtu)
