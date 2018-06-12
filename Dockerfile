FROM ubuntu:16.04
ADD . /pipeline/fusion/
RUN apt-get update && \
	apt-get upgrade -y && \
	apt-get install python -y && \
	apt-get install default-jre -y && \
	apt-get install python-numpy -y && \
	apt-get install python-pandas -y && \
	apt-get install python-sklearn -y && \
	apt-get clean && \
	cd /pipeline/fusion/pegasus/learn && \
	python train_model.py