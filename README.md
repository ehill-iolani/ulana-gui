# Ulana
***U***nicellular ***L***ong-read ***A***ssembly a***N***d ***A***nnotation

A bacterial genome assembly and annotation pipeline using Fast, HAC or SUP ONT basecalled data from MinION and Flongle flow cells. Updated with a GUI for ease of use.

ulana

1. vt. To plait, weave, knit, braid; plaiting, weaving. Also unala, nala, unana. Mea ulana ʻia, plaited or woven material, textile. Mea ulana lole, weaver (Isa. 38.12), loom. (PPN langa.)

# Installation

Building the docker image from the dockerfile:
```
cd path/to/directory/ulana-gui/
docker build -t ulana-gui:latest .
docker run --name=ulana-gui -p 8080:3838 -dt ulana-gui:latest
```

Using the dockerhub image:
```
docker pull ethill/ulana-gui:latest
docker run --name=ulana-gui -p 8080:3838 -dt ulana-gui:latest
```

Once the container is running navigate to http://localhost:8080 to use the application!
