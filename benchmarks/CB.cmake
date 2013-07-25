ADD_BENCHMARK ("CB" C/1D_isofrac/1d_isofrac_AS "OGS_FEM" 1
	C/1D_isofrac/1d_isofrac_AS_ply_OUT_LINE_t0.tec)

ADD_BENCHMARK ("CB" C/1D_isofrac/1d_isofrac "OGS_FEM" 1
	C/1D_isofrac/1d_isofrac_ply_OUT_LINE_t0.tec)

# Haibing is responsible for this benchmark.
#ADD_BENCHMARK ("CB" C/1d_xylene_degradation/h2_line "OGS_FEM" 1
#	C/1d_xylene_degradation/h2_line_domain_line.tec)

ADD_BENCHMARK ("CB" GROUNDWATER_FLOW/Transient_flow/trans_bd_homo "OGS_FEM" 1
	GROUNDWATER_FLOW/Transient_flow/trans_bd_homo_domain_GROUNDWATER_FLOW_quad.tec)

ADD_BENCHMARK ("CB_LONG" C/FG_3ports/rt1 "OGS_FEM" 1
	C/FG_3ports/rt1_domain_quad.tec)

ADD_BENCHMARK ("CB_LONG" C/hetk+n+restart/2D1P_transport "OGS_FEM" 1
	C/hetk+n+restart/2D1P_transport_domain_quad.tec
	C/hetk+n+restart/2D1P_transport_time_POINT10.tec
	C/hetk+n+restart/2D1P_transport_time_POINT11.tec)
