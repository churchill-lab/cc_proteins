FROM python:3.6-slim
LABEL maintainer="Matthew Vincent <mattjvincent@gmail.com>" \
	  version="1.0"

RUN apt-get update && \
    apt-get -y install gcc

ENV INSTALL_PATH /app/cc_proteins
RUN mkdir -p $INSTALL_PATH

WORKDIR $INSTALL_PATH

COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

COPY . .
RUN pip install --editable .

# install clustalo
RUN cp /app/cc_proteins/ext/clustalo-1.2.4-Ubuntu-x86_64 /usr/bin/clustalo-1.2.4 && \
    ln -s /usr/bin/clustalo-1.2.4 /usr/bin/clustalo

CMD gunicorn -c "python:config.gunicorn" --reload "mmc.app:create_app()"
