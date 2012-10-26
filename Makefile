wf.png: ho_gswv
	./ho_gswv > wf.txt
	python plot.py

clean:
	rm wf.png wf.txt ho_gswv
