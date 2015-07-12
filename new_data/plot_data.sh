#!/bin/bash

for i in $(ls -d */); do

    D_dir=${i}
    D_tmp=${D_dir#*=}
    D=${D_tmp%/*}

    cd ${i}

    for j in $(ls -d */); do

	Bo_dir=${j}
	Bo_tmp=${Bo_dir#*=}
	Bo=${Bo_tmp%/*}

	prod=$(echo $D*$Bo | bc -l)

	if [[ $prod == 1* ]]; 
	then
	    t_step=0.01
	else
	    temp=0.01
	    t_step=$(echo $temp*$prod | bc -l)
	fi
	cd ${j}

	for k in $(ls -d */); do

	    cd ${k}

	    rm *.png

	    gnuplot -e "traj_plot='sphere.eps'" -e "traj_data='sphere.dat'" ../../../plot_data.p

	    for f in interf_config*; do 
		filename=${f}
		tmp=${filename#*g}
		num=${tmp%.*}
		size=${#num}

		time=$(echo $t_step*$num | bc -l)

		if [ $size -eq "1" ]; then
		    plot_file="config0000"${num}
		elif [ $size -eq "2" ]; then
		    plot_file="config000"${num}
		elif [ $size -eq "3" ]; then
		    plot_file="config00"${num}
		elif [ $size -eq "4" ]; then
		    plot_file="config0"${num}
		elif [ $size -eq "5" ]; then 
		    plot_file="config"${num}
		fi

		sphere_data="sphere_config"${num}".dat"

		gnuplot -e "config_plot='${plot_file}.eps'" -e "interf_dat='${f}'" -e "sphere_dat='${sphere_data}'" -e "time='${time}'" ../../../plot_config.p

		convert -density 300 ${plot_file}".eps" -flatten ${plot_file}".png"

		convert ${plot_file}".png" -rotate 90 ${plot_file}".png"

		convert -delay 100 -loop 0 *.png config.gif

		rm ${plot_file}".eps"

	    done

	    convert -density 300 "sphere.eps" -flatten "sphere.png"
	    
	    convert "sphere.png" -rotate 90 "sphere.png"

	    echo $i, $j $k
	    cd ..
	done

	cd ..

    done

    cd ..
   
done
