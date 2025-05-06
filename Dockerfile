FROM --platform=$BUILDPLATFORM ubuntu:noble

RUN apt-get update && apt-get -y upgrade

RUN apt-get install -y build-essential cmake git time \
       ssh g++ r-cran-devtools samtools neofetch locales

RUN apt --purge autoremove && apt-get clean

RUN localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG=en_US.utf8

RUN R -e 'install.packages("ggmuller")'
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'install.packages("patchwork")'
RUN R -e 'install.packages("optparse")'
RUN R -e 'install.packages("ggrepel")'
RUN R -e 'install.packages("ggpubr")'
RUN R -e 'install.packages("bench")'
RUN R -e 'devtools::install_github("caravagnalab/ProCESS@1.0.1")'

ENV DEBIAN_FRONTEND=noninteractive
ARG USER_ID=ProCESS
RUN useradd -m $USER_ID
USER $USER_ID
WORKDIR /home/$USER_ID
ENV HOME=/home/$USER_ID

CMD ["neofetch"]
CMD ["time", "--version"]
