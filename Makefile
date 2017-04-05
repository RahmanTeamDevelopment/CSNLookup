
setup:
	echo "Setting up CSNLookup v1.0.0"
	virtualenv env
	source env/bin/activate; pip install -r requirements.txt
	mkdir submissions
	mkdir transdbs
	tar -zxvf CAVA-1.2.0.tar.gz
	mv CAVA-1.2.0 cava
	source env/bin/activate; cd cava; ./install.sh
	echo "Finished setting up CSNLookup v1.0.0"

deploy:
	eb create csn-lookup-env --timeout=5000 --verbose --debug
