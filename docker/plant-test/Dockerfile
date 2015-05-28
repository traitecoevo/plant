FROM r-base

RUN apt-get update \
  && apt-get install -y --no-install-recommends  \
    ca-certificates \
    curl \
    git \
    libcurl4-openssl-dev \
    libicu-dev \
    ssh \
  && apt-get clean

RUN install2.r --error  \
    BB \
    BH \
    deSolve \
    devtools \
    nleqslv \
    numDeriv \
    R6 \
    Rcpp \
    testthat

RUN installGithub.r  \
    traitecoevo/dockertest \
    richfitz/RcppR6

RUN r -e 'dockertest:::copy_scripts_dir("/usr/local/bin")'

RUN mkdir /root/plant \
  && echo "clone.sh /root/plant" >> /root/.bashrc \
  && echo "system('clone.sh /root/plant')" > /root/.Rprofile \
  && echo "system('clone.sh /root/plant')" > /root/.littler.r

WORKDIR /root/plant

CMD ["bash"]
