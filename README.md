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

# Citation 
You can find the Ulana paper [here](https://www.liebertpub.com/doi/10.1089/ast.2023.0072)

If you use this pipeline please use this citation:

Prescott, R. D., Chan, Y. L., Tong, E. J., Bunn, F., Onouye, C. T., Handel, C., Lo, C.-C., Davenport, K., Johnson, S., Flynn, M., Saito, J. A., Lee, H., Wong, K., Lawson, B. N., Hiura, K., Sager, K., Sadones, M., Hill, E. C., Esibill, D., … Donachie, S. P. (2023a). Bridging Place-based astrobiology education with genomics, including descriptions of three novel bacterial species isolated from Mars analog sites of cultural relevance. Astrobiology, 23(12), 1348–1367. https://doi.org/10.1089/ast.2023.0072 
