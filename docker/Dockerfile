FROM drone/rlang

WORKDIR /home/ubuntu
USER ubuntu

# We need some system software:
RUN sudo apt-get -y install libgsl0-dev

# And install some deps too.  Ideally we could just do:
#   RUN sudo travis-tool.sh install_deps
# but that requires getting the DESCRIPTION file into /home/ubuntu
#
# So for now, just doing it manually:
RUN sudo travis-tool.sh install_r_binary numDeriv deSolve digest
