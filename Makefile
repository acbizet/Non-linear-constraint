.PHONY: install exact50 analyses clean

install:
	python -m pip install -e .

exact50:
	python scripts/01_generate_exact_series.py --qmax 50 --progress

analyses:
	python scripts/reproduce_analysis.py

clean:
	rm -rf results/generated
