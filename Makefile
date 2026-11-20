.PHONY: run ci test clean

run:
	Rscript run_pipeline.R --config config/config.yaml

ci:
	Rscript run_pipeline.R --config config/ci.yaml

test:
	Rscript -e "testthat::test_dir('tests/testthat', stop_on_failure = TRUE)"

clean:
	rm -rf results/*.tsv results/*.rds figures/*.png data/raw/*.gz data/raw/*.txt data/raw/*.tsv data/raw/*.json
