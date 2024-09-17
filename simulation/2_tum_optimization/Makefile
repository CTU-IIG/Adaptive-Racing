all: Ruudskogen CG_Speedway_1


Ruudskogen: Ruudskogen_MCP Ruudskogen_MCP_IQP

Ruudskogen_MCP: outputs/ruuds_mincurv_f110.csv

Ruudskogen_MCP_IQP: outputs/ruuds_mincurv_iqp_f110.csv

outputs/ruuds_mincurv_f110.csv:
	./run_optimization.py3 --model params/f110.ini --opt_type mincurv --inputs_folder ../global_racetrajectory_optimization --output ruuds_mincurv_f110 --save-plot ../0_map_generation/Ruudskogen

outputs/ruuds_mincurv_iqp_f110.csv:
	./run_optimization.py3 --model params/f110.ini --opt_type mincurv_iqp --inputs_folder ../global_racetrajectory_optimization --output ruuds_mincurv_iqp_f110 --save-plot ../0_map_generation/Ruudskogen


CG_Speedway_1: CG_Speedway_1_MCP CG_Speedway_1_MCP_IQP

CG_Speedway_1_MCP: outputs/cgs_mincurv_f110.csv

CG_Speedway_1_MCP_IQP: outputs/cgs_mincurv_iqp_f110.csv

outputs/cgs_mincurv_f110.csv:
	./run_optimization.py3 --model params/f110.ini --opt_type mincurv --inputs_folder ../global_racetrajectory_optimization --output cgs_mincurv_f110 --save-plot ../0_map_generation/CG_Speedway_1

outputs/cgs_mincurv_iqp_f110.csv:
	./run_optimization.py3 --model params/f110.ini --opt_type mincurv_iqp --inputs_folder ../global_racetrajectory_optimization --output cgs_mincurv_iqp_f110 --save-plot ../0_map_generation/CG_Speedway_1
