# CSN_lab5: Evaluating Community Detection Algorithms
The objective of this project is to explore a variety of community detection algorithms, and see some of the _strengths_ and _weaknesses_ for all of them. 
In order to do so, we compute 4 metrics - triangle participation ratio, modularity, conductance, expansion - over the obtained communities on 4 different graphs. 
Two of them are provided by the __igraphdata__ package, while the other two where artificially created by merging networks in which members are tightly connected, 
and communities with less within-group connections.
- _Macaque Network_: Visuotactile brain areas and connections of the macaque monkey.
- _Hospital Network_: Records of contacts among patients and various types of health care workers in
the geriatric unit of a hospital in Lyon, France, in 2010.
- _Artificial Network from ER Graphs_.
- _Artificial Network from Complete Graphs_.
Finally, to conclude, we have applied those concept on a bigger network Wikipedia.gml.
