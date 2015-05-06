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

RUN mkdir /root/tree2 \
  && echo "clone.sh /root/tree2" >> /root/.bashrc \
  && echo "system('clone.sh /root/tree2')" > /root/.Rprofile \
  && echo "system('clone.sh /root/tree2')" > /root/.littler.r

WORKDIR /root/tree2

CMD ["bash"]
