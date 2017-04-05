
setup:
	echo "Setting up CSNLookup v1.0.0"
	mkdir submissions
	mkdir transdbs
	tar -zxvf CAVA-1.2.0.tar.gz
	mv CAVA-1.2.0 cava
	cd cava; ./install.sh
	echo "Finished setting up CSNLookup v1.0.0"

deploy:
	eb create csn-lookup-env --timeout=5000 --verbose --debug
