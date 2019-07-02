FROM r-base

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    libcurl4-openssl-dev \
    libicu-dev \
    libssh2-1-dev \
    libssl-dev \
    libxml2-dev \
    python-dev \
    ssh \
    zlib1g-dev \
  && apt-get clean

RUN install2.r --error  \
    BB \
    BH \
    crayon \
    deldir \
    deSolve \
    devtools \
    pkgbuild \
    knitr \
    nleqslv \
    numDeriv \
    R6 \
    Rcpp \
    rPython

RUN installGithub.r  \
#    traitecoevo/callr \
#    traitecoevo/dockertest \
    richfitz/RcppR6 \
    hadley/testthat \
#    richfitz/grader \
    smbache/loggr \
    jimhester/covr

RUN r -e 'dockertest:::copy_scripts_dir("/usr/local/bin")'

COPY local /local

RUN R CMD INSTALL /local/plant.ml \
  && rm -rf /local

RUN mkdir /root/plant \
  && echo "clone.sh /root/plant" >> /root/.bashrc \
  && echo "system('clone.sh /root/plant')" > /root/.Rprofile \
  && echo "system('clone.sh /root/plant')" > /root/.littler.r

WORKDIR /root/plant

CMD ["bash"]
