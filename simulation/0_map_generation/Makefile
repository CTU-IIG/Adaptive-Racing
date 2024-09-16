# (Example) Generate map
all: Ruudskogen CG_Speedway_1

Ruudskogen: Ruudskogen.csv ruuds_start_points.npy

CG_Speedway_1: CG_Speedway_1.csv cgs_start_points.npy

Ruudskogen.png:
	convert Ruudskogen.svg -scale 7351x4395 Ruudskogen.png

Ruudskogen.pgm:
	convert Ruudskogen.svg -scale 735x439 -depth 8 Ruudskogen.pgm

ruuds_start_points.npy: Ruudskogen.png
	ng_generate_data -f Ruudskogen.png --resolution 0.042 --resize 0.1 --scale 0.1 --ox 0 --oy 0 --prefix ruuds_

Ruudskogen.csv: Ruudskogen.pgm
	python3 map_converter.py3 Ruudskogen

CG_Speedway_1.png:
	convert CG_Speedway_1.svg -scale 1505x1425 CG_Speedway_1.png

CG_Speedway_1.pgm:
	convert CG_Speedway_1.svg -scale 301x285 -depth 8 CG_Speedway_1.pgm

cgs_start_points.npy: CG_Speedway_1.png
	ng_generate_data -f CG_Speedway_1.png --resolution 0.25862069 --resize 0.2 --scale 0.1 --ox 0 --oy 0 --prefix cgs_

CG_Speedway_1.csv: CG_Speedway_1.pgm
	python3 map_converter.py3 CG_Speedway_1

clean:
	rm -f *.png

clean-all: clean
	rm -f *.pgm *.npy *.csv
