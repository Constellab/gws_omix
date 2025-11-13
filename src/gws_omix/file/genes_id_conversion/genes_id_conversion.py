# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import shlex
from pathlib import Path
from typing import Final
import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import sys

from gws_core import(
    ConfigParams, ConfigSpecs, File, FloatParam, InputSpec, InputSpecs,
    OutputSpec, OutputSpecs, PlotlyResource, ResourceSet, ShellProxy,
    StrParam, TableImporter, Task, TaskInputs, TaskOutputs, task_decorator,TableImporter,Table)

from .genes_id_conversion_env import GenesidConversionShellProxyHelper

SCIENTIFIC_NAMES = [
"Daucus carota subsp. sativus","Bombus terrestris","Teladorsagia circumcincta","Ananas comosus",
"Triticum aestivum","Juglans regia","Angiostrongylus cantonensis","Bactrocera latifrons","Drosophila ficusphila",
"Petromyzon marinus","Trichinella murrelli","Bicyclus anynana","Pyrenophora teres f. teres 0-1","Linhomoeus sp. GSCO2_2",
"Poecilia latipinna","Trichobilharzia regenti","Cajanus cajan","Ciona intestinalis","Hermetia illucens","Schistosoma intercalatum",
"Caenorhabditis nigoni","Steinernema carpocapsae","Atriophallophorus winterbourni","Mesocricetus auratus","Nematostella vectensis",
"Parastrongyloides trichosuri","Erpetoichthys calabaricus","Anneissia japonica","Vitis vinifera","Schistosoma curassoni","Jaculus jaculus",
"Anopheles atroparvus","Anopheles dirus","Parapristionchus giblindavisi","Schistosoma haematobium","Paragonimus kellicotti",
"Fusarium graminearum","Laticauda laticaudata","Thelazia callipaeda","Zymoseptoria tritici IPO323","Schistosoma rodhaini",
"Phytophthora lateralis MPF4","Arabis alpina","Macaca nemestrina","Lupinus angustifolius","Solenopsis invicta","Dimorphilus gyrociliatus",
"Lutzomyia longipalpis","Elaeophora elaphi","Bemisia tabaci","Macrosteles quadrilineatus","Marmota marmota marmota","Ailuropoda melanoleuca",
"Microctonus aethiopoides","Oryzias melastigma","Globodera rostochiensis","Mizuhopecten yessoensis","Echinococcus granulosus",
"Varroa destructor","Oryza glaberrima","Oryzias latipes","Oryza sativa aus subgroup","Schistosoma margrebowiei","Drosophila innubila",
"Plasmodium vivax","Manduca sexta","Heterodera schachtii","Vicia faba","Diploscapter pachys","Lingula anatina","Panthera pardus",
"Aspergillus terreus NIH2624","Dreissena polymorpha","Vigna radiata var. radiata","Sus scrofa","Caenorhabditis quiockensis",
"Anas platyrhynchos","Caenorhabditis waitukubuli","Anoplophora glabripennis","Cynoglossus semilaevis","Aphelenchoides fujianensis",
"Steinernema feltiae","Cylicocyclus nassatus","Drosophila hydei","Vombatus ursinus","Takifugu rubripes","Strongyloides papillosus",
"Secale cereale","Felis catus","Dicentrarchus labrax","Naja naja","Machimus atricapillus","Oncorhynchus mykiss","Caenorhabditis tropicalis",
"Anolis carolinensis","Galdieria sulphuraria","Arabidopsis thaliana","Oryza longistaminata","Spirometra erinaceieuropaei",
"Levipalatum texanum","Paragonimus westermani","Trichuris suis","Vicugna pacos","Oryza sativa Indica Group","Struthio camelus australis",
"Lottia gigantea","Oreochromis niloticus","Ditylenchus destructor","Paragonimus heterotremus","Quercus lobata","Rhynchonema sp. JSB1_4",
"Eremothecium gossypii ATCC 10895","Drosophila mojavensis","Ditylenchus dipsaci","Heligmosomoides polygyrus","Lolium perenne",
"Aquila chrysaetos chrysaetos","Brassica napus","Actinia tenebrosa","Theobroma cacao","Pediculus humanus corporis",
"Danaus plexippus plexippus","Globodera pallida","Caenorhabditis zanzibari","Lucilia cuprina","Fasciola hepatica","Xenopus tropicalis",
"Fusarium fujikuroi","Hydra vulgaris","Carassius auratus","Schistocerca serialis cubense","Microtus ochrogaster","Stegodyphus mimosarum",
"Glossina brevipalpis","Drosophila persimilis","Pan troglodytes","Beauveria bassiana","Thalassiosira pseudonana CCMP1335",
"Drosophila guanche","Melitaea cinxia","Trichinella pseudospiralis","Caenorhabditis briggsae","Drosophila yakuba","Brachypodium distachyon",
"Heterocephalus glaber","Acyrthosiphon pisum","Caenorhabditis elegans","Meloidogyne incognita","Bursaphelenchus xylophilus",
"Trichinella patagoniensis","Paramecium tetraurelia","Setaria viridis","Helianthus annuus","Echinostoma caproni","Amphimedon queenslandica",
"Lytechinus pictus","Corylus avellana","Daktulosphaira vitifoliae","Anthonomus grandis grandis","Dothistroma septosporum NZE10",
"Panagrolaimus sp. JU765","Terrapene triunguis","Anopheles arabiensis","Sarcoptes scabiei","Caenorhabditis tribulationis",
"Pollicipes pollicipes","Cherax quadricarinatus","Dioscorea cayenensis subsp. rotundata","Schistosoma mattheei","Amborella trichopoda",
"Neogale vison","Strongylocentrotus purpuratus","Canis lupus familiaris","Haliotis rufescens","Adelges cooleyi","Centruroides sculpturatus",
"Onchocerca ochengi","Oncorhynchus tshawytscha","Oryza sativa Japonica Group","Schistocerca gregaria","Varroa jacobsoni",
"Myopa tessellatipennis","Eriocheir sinensis","Phaseolus vulgaris","Litomosoides sigmodontis","Aspergillus nidulans FGSC A4",
"Syphacia muris","Drosophila virilis","Phanerochaete chrysosporium RP-78","Ixodes scapularis","Bombyx mandarina","Taenia multiceps",
"Strongylus vulgaris","Salmo trutta","Glossina morsitans morsitans","Anopheles melas","Tetranychus urticae","Aphidius gifuensis",
"Glossina palpalis gambiensis","Rhodnius prolixus","Drosophila simulans","Cavia porcellus","Acropora millepora","Eragrostis tef",
"Necator americanus","Schistosoma bovis","Pristionchus arcanus","Schistosoma turkestanicum","Mus spretus","Panonychus citri",
"Physcomitrium patens","Priapulus caudatus","Diplogasteroides magnus","Digitaria exilis","Trichinella nativa","Drosophila eugracilis",
"Catagonus wagneri","Angiostrongylus vasorum","Cercopithifilaria johnstoni","Mercenaria mercenaria","Cynara cardunculus var. scolymus",
"Cryptococcus deneoformans JEC21","Tetrahymena thermophila SB210","Halicephalobus mephisto","Cannabis sativa","Rhinolophus ferrumequinum",
"Drosophila sechellia","Athalia rosae","Pristionchus mayeri","Auanema sp. JU1783","Otolemur garnettii","Verticillium dahliae VdLs.17",
"Panagrolaimus sp. ES5","Solanum tuberosum","Drosophila arizonae","Phocoena sinus","Enoplolaimus lenunculus","Pristionchus maxplancki",
"Nakaseomyces glabratus CBS 138","Physeter macrocephalus","Aethina tumida","Avena sativa","Caenorhabditis inopinata",
"Anas platyrhynchos platyrhynchos","Caenorhabditis uteleia","Trichuris muris","Hippocampus comes","Salvator merianae",
"Diuraphis noxia","Canis lupus dingo","Trichogramma pretiosum","Drosophila biarmipes","Galleria mellonella","Cyprinus carpio carpio",
"Acanthaster planci","Macrostomum lignano","Phytophthora sojae","Trissonchulus sp. WLG1_4","Parastagonospora nodorum SN15",
"Saccharomyces cerevisiae S288C","Schistosoma japonicum","Seriola dumerili","Tursiops truncatus","Strongyloides ratti","Cyprinodon variegatus",
"Helicoverpa zea","Hordeum vulgare subsp. vulgare","Haplochromis burtoni","Dictyocaulus viviparus","Hyalella azteca","Macaca fascicularis",
"Coremacera marginata","Limnephilus lunatus","Marchantia polymorpha","Entamoeba histolytica HM-1:IMSS","Tetraodon nigroviridis",
"Lepisosteus oculatus","Lepeophtheirus salmonis","Dendronephthya gigantea","Mayetiola destructor","Plasmodium chabaudi chabaudi",
"Trichinella nelsoni","Teleopsis dalmanni","Gossypium raimondii","Loa loa","Nasonia vitripennis","Oryza sativa tropical japonica subgroup",
"Caenorhabditis parvicauda","Taenia asiatica","Diabrotica virgifera virgifera","Eutrema salsugineum","Larimichthys crocea",
"Rhopalosiphum maidis","Candidozyma duobushaemuli","Schistocerca nitens","Sorghum bicolor","Heterobilharzia americana","Apis mellifera",
"Owenia fusiformis","Agrilus planipennis","Beta vulgaris subsp. vulgaris","Cebus imitator","Meloidogyne graminicola","Ascaris lumbricoides",
"Prunus dulcis","Pomacea canaliculata","Meloidogyne enterolobii","Halicephalobus sp. NKZ332","Halyomorpha halys","Albugo laibachii Nc14","Drosophila busckii","Adineta vaga","Saccharum spontaneum","Xiphophorus couchianus","Opisthorchis felineus","Meloidogyne arenaria","Lactuca sativa","Aspergillus flavus NRRL3357","Glossina pallidipes","Chelonus insularis","Chondrus crispus","Saccoglossus kowalevskii","Sorex araneus","Pyricularia oryzae 70-15","Candidozyma haemuli","Anopheles funestus","Ichneumon xanthorius","Rhagoletis pomonella","Drosophila pseudoobscura","Trialeurodes vaporariorum","Aphelenchoides besseyi","Amphiprion ocellaris","Glyphotaelius pellucidus","Ixodes persulcatus","Caenorhabditis panamensis","Salmo salar","Leishmania major strain Friedlin","Actinia equina","Panagrellus redivivus","Steinernema glaseri","Manihot esculenta","Leersia perrieri","Sciurus vulgaris","Nymphaea colorata","Heterodera glycines","Drosophila elegans","Toxocara canis","Candida parapsilosis CDC317","Solanum lycopersicum","Bactrocera dorsalis","Prunus avium","Trichinella sp. T9","Limnephilus marmoratus","Wuchereria bancrofti","Bison bison bison","Aspergillus fischeri NRRL 181","Micoletzkya japonica","Puccinia striiformis f. sp. tritici PST-130","Trypanosoma brucei","Erinaceus europaeus","Plasmodium falciparum 3D7","Trichinella spiralis","Bos mutus","Aspergillus clavatus NRRL 1","Clonorchis sinensis","Lineus longissimus","Anopheles christyi","Oryza glumipatula","Rhabditophanes sp. KR3021","Denticeps clupeoides","Galendromus occidentalis","Capsicum annuum","Phlebotomus papatasi","Schmidtea mediterranea","Geospiza fortis","Brugia malayi","Astyanax mexicanus","Dinothrombium tinctorium","Selaginella moellendorffii","Hyaloperonospora arabidopsidis Emoy2","Gigantopelta aegis","Globisporangium ultimum DAOM BR144","Pseudo-nitzschia multistriata","Ancylostoma ceylanicum","Myripristis murdjan","Serinus canaria","Gadus morhua","Citrullus lanatus","Opisthorchis viverrini","Ictidomys tridecemlineatus","Cyclopterus lumpus","Ctenocephalides felis","Electrophorus electricus","Clupea harengus","Populus trichocarpa","Plasmodium knowlesi strain H","Limnephilus rhombicus","Guillardia theta CCMP2712","Heliconius melpomene","Dibothriocephalus latus","Hucho hucho","Acrobeloides nanus","Mustela putorius furo","Nomascus leucogenys","Ustilago maydis","Prunus persica","Pectinophora gossypiella","Oryzias sinensis","Ooceraea biroi","Caenorhabditis remanei","Paragonimus skrjabini miyazakii","Meloidogyne floridensis","Harpegnathos saltator","Ascaris suum","Ovis aries","Orbicella faveolata","Globisporangium irregulare DAOM BR486","Bos indicus x Bos taurus","Panicum hallii","Pristionchus entomophagus","Caenorhabditis auriculariae","Phytopythium vexans DAOM BR484","Dendroctonus ponderosae","Gasterosteus aculeatus aculeatus","Schizosaccharomyces japonicus yFS275","Patiria miniata","Meloidogyne chitwoodi","Biomphalaria glabrata","Oryza rufipogon","Aedes albopictus","Dracunculus medinensis","Gyrodactylus bullatarudis","Medicago truncatula","Danio rerio","Carlito syrichta","Sabatieria punctata","Saimiri boliviensis boliviensis","Octopus bimaculoides","Panagrolaimus davidi","Microctonus hyperodae","Musa acuminata subsp. malaccensis","Leptotrombidium deliense","Bos grunniens","Brassica juncea","Sesamum indicum","Anopheles quadriannulatus","Oesophagostomum dentatum","Onchocerca flexuosa","Ficedula albicollis","Lytechinus variegatus","Equus asinus","Aspergillus fumigatus A1163","Globisporangium iwayamae DAOM BR242034","Ancylostoma duodenale","Phytophthora kernoviae 00238/432","Paramacrobiotus metropolitanus","Equus caballus","Microbotryum lychnidis-dioicae p1A1 Lamole","Zea mays","Drosophila albomicans","Phytophthora infestans T30-4","Drosophila takahashii","Phytophthora ramorum","Botrytis cinerea B05.10","Haemonchus contortus","Gyrodactylus salaris","Notechis scutatus","Steinernema hermaphroditum","Nannospalax galili","Crocodylus porosus","Dermatophagoides pteronyssinus","Brugia pahangi","Coturnix japonica","Ancylostoma caninum","Ursus maritimus","Panicum hallii var. hallii","Schizosaccharomyces octosporus yFS286","Penaeus vannamei","Bradynema listronoti","Drosophila navojoa","Parus major","Dictyostelium discoideum AX4","Macaca mulatta","Nicotiana attenuata","Arabidopsis lyrata subsp. lyrata","Heterorhabditis bacteriophora","Quercus suber","Exaiptasia diaphana","Stylophora pistillata","Fusarium verticillioides 7600","Dipodomys ordii","Portunus trituberculatus","Thrips palmi","Melampsora larici-populina 98AG31","Trichinella sp. T8","Habropoda laboriosa","Camelina sativa","Ostrea edulis","Pteropus vampyrus","Anopheles gambiae","Malus domestica","Bombus huntii","Tribolium madens","Drosophila grimshawi","Cottoperca gobio","Theristus sp. LFF4_11","Leptinotarsa decemlineata","Cimex lectularius","Camponotus floridanus","Anopheles sinensis","Rodentolepis nana","Melanaphis sacchari","Scophthalmus maximus","Chelonoidis abingdonii","Pyricularia oryzae","Paralinhomoeus sp. GSCO2_6","Callorhinchus milii","Kalanchoe fedtschenkoi","Trichoderma virens Gv29-8","Helobdella robusta","Drosophila rhopaloa","Puccinia graminis f. sp. tritici 04KEN156/4","Fusarium oxysporum f. sp. lycopersici 4287","Papio anubis","Myotis lucifugus","Setaria italica","Dicrocoelium dendriticum","Olea europaea subsp. europaea","Gaeumannomyces tritici R3-111a-1","Fusarium pseudograminearum CS3096","Tupaia belangeri","Sitophilus oryzae","Pythium arrhenomanes ATCC 12531","Pocillopora damicornis","Haemaphysalis longicornis","Apis dorsata","Haliotis rubra","Microlaimidae sp. YZB2_3","Strigops habroptila","Stegodyphus dumicola","Protopolystoma xenopodis","Anopheles culicifacies","Dasypus novemcinctus","Maylandia zebra","Mesodorylaimus sp. YZB2_4","Strongyloides venezuelensis","Brugia timori","Gopherus evgoodei","Candidozyma pseudohaemuli","Brassica rapa","Triticum dicoccoides","Glycine soja","Acanthochromis polyacanthus","Pogonomyrmex barbatus","Parascaris univalens","Pristionchus pacificus","Bursaphelenchus okinawaensis","Mastacembelus armatus","Arabidopsis halleri subsp. gemmifera","Podarcis muralis","Amblyteles armatorius","Plectus sambesii","Gallus gallus","Mus caroli","Eucalyptus grandis","Procambarus clarkii","Drosophila gunungcola","Notamacropus eugenii","Anser brachyrhynchus","Zootermopsis nevadensis","Prolemur simus","Echinococcus multilocularis","Amphiprion percula","Octopus sinensis","Drosophila willistoni","Cercocebus atys","Clytia hemisphaerica","Dirofilaria immitis","Aspergillus oryzae RIB40","Rhinopithecus roxellana","Aphanomyces invadans","Poecilia formosa","Drosophila kikkawai","Homo sapiens","Homalodisca vitripennis","Fusarium culmorum CS7071","Aphanomyces astaci","Panagrolaimus superbus","Pelodiscus sinensis","Blumeria hordei DH14","Betta splendens","Daphnia pulicaria","Sander lucioperca","Glossina fuscipes","Anopheles farauti","Asparagus officinalis","Cricetulus griseus","Acanthocheilonema viteae","Candidozyma auris","Penaeus japonicus","Drosophila miranda","Brassica oleracea var. oleracea","Propithecus coquereli","Hyalomma asiaticum","Urocitellus parryii","Aspergillus fumigatus Af293","Giardia lamblia ATCC 50803","Anopheles maculatus","Triticum turgidum subsp. durum","Chara braunii","Cryptotermes secundus","Parascaris equorum","Coffea canephora","Aphelenchoides bicaudatus","Ficus carica","Vulpes vulpes","Anopheles stephensi","Anisakis simplex","Rhipicephalus microplus","Oryza sativa aromatic subgroup","Histoplasma ohiense (nom. inval.)","Hypsibius exemplaris","Verticillium dahliae JR2","Puccinia triticina 1-1 BBBD Race 1","Phascolarctos cinereus","Pecten maximus","Oppia nitens","Esox lucius","Bactrocera neohumeralis","Panagrolaimus sp. PS1159","Dermacentor silvarum","Nippostrongylus brasiliensis","Magnaporthiopsis poae ATCC 64411","Leptobrachium leishanense","Moschus moschiferus","Pomphorhynchus laevis","Uloborus diversus","Cervus hanglu yarkandensis","Ursus americanus","Drosophila pseudoobscura pseudoobscura","Aedes aegypti","Allodiplogaster sudhausi","Amphibalanus amphitrite","Brassica rapa subsp. trilocularis","Rosa chinensis","Aegilops umbellulata","Emiliania huxleyi CCMP1516","Trichoplax adhaerens","Anopheles coluzzii","Olea europaea var. sylvestris","Ciona savignyi","Polistes dominula","Stegastes partitus","Mnemiopsis leidyi","Oryctolagus cuniculus","Taenia saginata","Komagataella phaffii GS115","Corchorus capsularis","Tribolium castaneum","Anopheles merus","Sclerotinia sclerotiorum 1980 UF-70","Stomoxys calcitrans","Kryptolebias marmoratus","Caenorhabditis sinica","Callithrix jacchus","Plasmodium berghei ANKA","Daphnia carinata","Colletotrichum fructicola Nara gc5","Trichinella sp. T6","Caenorhabditis brenneri","Neodiprion lecontei","Panthera leo","Caenorhabditis latens","Cataglyphis hispanica","Choloepus hoffmanni","Labrus bergylta","Plenodomus lingam JN3","Corymbia citriodora subsp. variegata","Epsilonema sp. ZAB3_2","Nilaparvata lugens","Sipha flava","Papaver somniferum","Sparus aurata","Oryza brachyantha","Procavia capensis","Pristionchus fissidentatus","Drosophila melanogaster","Monodelphis domestica","Ipomoea triloba","Patella pellucida","Pygocentrus nattereri","Meloidogyne hapla","Triticum urartu","Ornithorhynchus anatinus","Toxoplasma gondii ME49","Microcebus murinus","Meloidogyne javanica","Aplysia californica","Glycine max","Chinchilla lanigera","Phaeodactylum tricornutum CCAP 1055/1","Penaeus monodon","Drosophila erecta","Sarcophilus harrisii","Caenorhabditis angaria","Gongylonema pulchrum","Polistes canadensis","Pristionchus japonicus","Schistocerca americana","Hymenolepis diminuta","Ochotona princeps","Angiostrongylus costaricensis","Schizosaccharomyces cryophilus OY26","Haemonchus placei","Atta cephalotes","Rhinopithecus bieti","Delphinapterus leucas","Monomorium pharaonis","Amyelois transitella","Candida albicans SC5314","Oryza barthii","Belgica antarctica","Mus musculus","Taeniopygia guttata","Trichoderma reesei QM6a","Magallana gigas","Balaenoptera musculus","Aotus nancymaae","Mesorhabditis belari","Diaphorina citri","Lates calcarifer","Nothobranchius furzeri","Drosophila subobscura","Daphnia pulex","Hymenolepis microstoma","Pistacia vera","Oryza punctata","Cyanidioschyzon merolae strain 10D","Drosophila suzukii","Mandrillus leucophaeus","Drosophila ananassae","Trichinella britovi","Bombus impatiens","Colletotrichum graminicola M1.001","Echinochloa crus-galli","Steinernema monticolum","Anopheles albimanus","Anabas testudineus","Trileptium ribeirensis","Colletotrichum orbiculare MAFF 240422","Trifolium pratense","Capra hircus","Capitella teleta","Oscheius tipulae","Diploscapter coronatus","Loxodonta africana","Crassostrea virginica","Panthera tigris altaica","Trichuris trichiura","Xiphophorus maculatus","Hofstenia miamia","Octodon degus","Schistosoma mansoni","Schistocerca cancellata","Zerene cesonia","Megaselia scalaris","Venturia canescens","Meleagris gallopavo","Eragrostis curvula","Camelus dromedarius","Drosophila subpulchrella","Bactrocera tryoni","Bunonema sp. RGD898","Schistosoma guineensis","Monodon monoceros","Actinidia chinensis var. chinensis","Megachile rotundata","Rhipicephalus sanguineus","Caenorhabditis sulstoni","Candida tropicalis MYA-3404","Citrus x clementina","Orchesella cincta","Mus spicilegus","Limulus polyphemus","Limnoperna fortunei","Leguminivora glycinivorella","Sphenodon punctatus","Folsomia candida","Drosophila santomea","Cucumis sativus","Colletotrichum higginsianum","Strigamia maritima","Scleropages formosus","Neodiprion pinetum","Culex quinquefasciatus","Chrysemys picta bellii","Asterias rubens","Oryza meridionalis","Ceratitis capitata","Pongo abelii","Neurospora crassa OR74A","Orussus abietinus","Drosophila mauritiana","Ptycholaimellus sp. GST1_10","Schistosoma spindalis","Strongyloides stercoralis","Drosophila teissieri","Trichobilharzia szidati","Enterobius vermicularis","Neolamprologus brichardi","Setaria digitata","Bombus vancouverensis nearcticus","Trichinella zimbabwensis","Cucumis melo","Apis florea","Hydractinia symbiolongicarpus","Thelohanellus kitauei","Schistocephalus solidus","Pseudonaja textilis","Pan paniscus","Chlorocebus sabaeus","Fundulus heteroclitus","Phytophthora nicotianae P1569","Branchiostoma lanceolatum","Mesocestoides corti","Pyrenophora tritici-repentis Pt-1C-BFP","Peromyscus maniculatus bairdii","Dufourea novaeangliae","Mesorhabditis spiculigera","Yarrowia lipolytica CLIB122","Fasciola gigantica","Vigna angularis","Trichinella papuae","Soboliphyme baturini","Mya arenaria","Cotesia glomerata","Aegilops tauschii subsp. strangulata","Linepithema humile","Seriola lalandi dorsalis","Homarus americanus","Schizosaccharomyces pombe 972h-","Amphilophus citrinellus","Oryzias javanicus","Dermacentor andersoni","Pisum sativum","Astatotilapia calliptera","Glossina austeni","Gorilla gorilla gorilla","Hydatigera taeniaeformis","Anopheles minimus","Fraxinus excelsior","Caenorhabditis japonica","Cylicostephanus goldi","Polistes fuscatus","Steinernema scapterisci","Anopheles darlingi","Helicoverpa armigera","Caenorhabditis becei","Penaeus chinensis","Pundamilia nyererei","Drosophila obscura","Puccinia graminis f. sp. tritici CRL 75-36-700-3","Latimeria chalumnae","Anopheles epiroticus","Phlebotomus perniciosus","Parasteatoda tepidariorum","Culicoides sonorensis","Fusarium vanettenii 77-13-4","Chlamydomonas reinhardtii","Homarus gammarus","Echinococcus canadensis","Taenia solium","Chenopodium quinoa","Oryza nivara","Patella vulgata","Parelaphostrongylus tenuis","Koerneria luziae","Mus pahari","Tuber melanosporum Mel28","Daphnia magna","Bos taurus","Pristionchus exspectatus","Bigelowiella natans CCMP2755","Eufriesea mexicana","Ostreococcus lucimarinus CCE9901","Copidosoma floridanum","Tigriopus californicus","Echinococcus oligarthrus","Schistocerca piceifrons","Eurytemora affinis","Triticum spelta","Echinops telfairi","Trissonchulus latispiculum","Aspergillus niger","Pythium aphanidermatum DAOM BR444","Onthophagus taurus","Ancistrocerus nigricornis","Bombyx mori","Oncorhynchus kisutch","Fasciolopsis buskii","Eptatretus burgeri","Drosophila bipectinata","Acromyrmex echinatior","Onchocerca volvulus","Vigna unguiculata","Rattus norvegicus","Romanomermis culicivorax","Caenorhabditis bovis","Sporisorium reilianum SRZ2"
]

TARGET_NAMES = ["AFFY_HC_G110","AFFY_HG_FOCUS","AFFY_HG_U133_PLUS_2","AFFY_HG_U133A_2","AFFY_HG_U133B","AFFY_HG_U95A","AFFY_HG_U95AV2" "AFFY_HG_U95B","AFFY_HG_U95C","AFFY_HG_U95D","AFFY_HG_U95E","AFFY_HT_HG_U133_PLUS_PM","AFFY_HTA_2_0","AFFY_HUGENEFL","AFFY_PRIMEVIEW","AFFY_U133_X3P","AGILENT_CGH_44B","AGILENT_GPL19072","AGILENT_GPL26966","AGILENT_GPL6848","AGILENT_SUREPRINT_G3_GE_8X60K","AGILENT_SUREPRINT_G3_GE_8X60K_V2","AGILENT_WHOLEGENOME","AGILENT_WHOLEGENOME_4X44K_V1","AGILENT_WHOLEGENOME_4X44K_V2","ARRAYEXPRESS","CCDS","CCDS_ACC","CHEMBL","CODELINK_CODELINK","DBASS3","DBASS5","EMBL","ENS_LRG_GENE","ENS_LRG_TRANSCRIPT","ENSG","ENSP","ENST","ENTREZGENE","ENTREZGENE_TRANS_NAME","GENECARDS","GO","GOSLIM_GOA","HGNC","HGNC_ACC","HGNC_TRANS_NAME","HPA","HPA_ACC","ILLUMINA_HUMANREF_8_V3","ILLUMINA_HUMANWG_6_V3","MEROPS","MIM_GENE","MIM_MORBID","MIRBASE","MIRBASE_ACC","MIRBASE_TRANS_NAME","PDB","PHALANX_ONEARRAY","PROTEIN_ID","PROTEIN_ID_ACC","REACTOME","REACTOME_GENE","REACTOME_TRANSCRIPT","REFSEQ_MRNA","REFSEQ_MRNA_ACC","REFSEQ_MRNA_PREDICTED","REFSEQ_MRNA_PREDICTED_ACC","REFSEQ_NCRNA","REFSEQ_NCRNA_ACC","REFSEQ_NCRNA_PREDICTED","REFSEQ_NCRNA_PREDICTED_ACC","REFSEQ_PEPTIDE","REFSEQ_PEPTIDE_ACC","REFSEQ_PEPTIDE_PREDICTED","REFSEQ_PEPTIDE_PREDICTED_ACC","RFAM","RFAM_ACC","RFAM_TRANS_NAME","RNACENTRAL","UCSC","UNIPARC","UNIPROT_GN","UNIPROT_GN_ACC","UNIPROT_ISOFORM","UNIPROTSPTREMBL","UNIPROTSPTREMBL_ACC","UNIPROTSWISSPROT","UNIPROTSWISSPROT_ACC","WIKIGENE","AFFY_HUEX_1_0_ST_V2","AFFY_HUGENE_1_0_ST_V1","AFFY_HUGENE_2_0_ST_V1","AFFY_HUGENE_2_1_ST_V1","BIOGRID","DBASS3_ACC","DBASS5_ACC","ENTREZGENE_ACC","GENECARDS_ACC","MIM_GENE_ACC","MIM_MORBID_ACC","WIKIGENE_ACC"]

NUMERIC_IDS_TREATED = [
    "auto",
    "AFFY_HUEX_1_0_ST_V2",
    "AFFY_HUGENE_1_0_ST_V1",
    "AFFY_HUGENE_2_0_ST_V1",
    "AFFY_HUGENE_2_1_ST_V1",
    "BIOGRID",
    "DBASS3_ACC",
    "DBASS5_ACC",
    "ENTREZGENE_ACC",
    "GENECARDS_ACC",
    "MIM_GENE_ACC",
    "MIM_MORBID_ACC",
    "WIKIGENE_ACC",
]


@task_decorator(
    "IDConvertWithGProfiler",
    human_name="Gene ID conversion",
    short_description="Convert gene IDs to a selected target namespace using g:Profiler.",
)
class IDConvertTask(Task):
    """
    This task (task: OmiX – Gene ID conversion with g:Profiler) performs identifier harmonization on a user-provided list of gene IDs. It uses
    the g:Profiler g:Convert service to translate input identifiers (e.g., gene symbols, Ensembl IDs, Entrez IDs, UniProt accessions) into a
    selected target namespace.
    See the full documentation on Cnostellab Community :https://constellab.community/bricks/gws_omix/latest/doc/use-cases/gene-id-conversion/
    """

    input_specs: Final[InputSpecs] = InputSpecs({
        "table_file": InputSpec(
            File,
            human_name="Input table (CSV/TSV)",
            short_description="Must contain the column to convert.",
        ),
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "converted_table": OutputSpec(Table, human_name="Converted IDs"),
        "annotated_file": OutputSpec(Table, human_name="Annotated input (+ g:Profiler columns)"),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "organism_name": StrParam(
            allowed_values=SCIENTIFIC_NAMES,
            short_description="Scientific name or g:Profiler code (e.g., 'Homo sapiens').",
        ),
        "id_column": StrParam(
            default_value="",
            short_description="Column to convert (leave blank to use the first column).",
        ),
        "target_namespace": StrParam(
            allowed_values=TARGET_NAMES,
            default_value="ENSG",
            short_description="g:Profiler target namespace.",
        ),
        "numeric_namespace": StrParam(
            allowed_values=NUMERIC_IDS_TREATED,
            default_value="auto",
            short_description="How to treat bare numeric IDs (e.g., ENTREZGENE_ACC). Recommended when IDs are numbers; otherwise set 'auto'.",
        ),
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_genes_id_conversion.py",
    )

    def run(self, p: "ConfigParams", ins: "TaskInputs") -> "TaskOutputs":
        infile      = ins["table_file"].path
        organism    = p["organism_name"].strip()
        id_column   = (p["id_column"] or "").strip()
        target_ns   = p["target_namespace"].strip()
        numeric_ns  = (p.get("numeric_namespace", "") or "").strip()

        shell: "ShellProxy" = GenesidConversionShellProxyHelper.create_proxy(self.message_dispatcher)
        work = shell.working_dir
        out_prefix = "ID_CONVERT"  # fixed prefix for outputs

        # Build argv safely, then quote as a single shell command for shell_mode=True
        argv = [
            "python3", self.python_file_path,
            "--infile", infile,
            "--organism", organism,
            "--target_ns", target_ns,
            "--out_prefix", out_prefix,
        ]
        if id_column:
            argv += ["--id_column", id_column]
        if numeric_ns:
            argv += ["--numeric_ns", numeric_ns]

        cmd = " ".join(shlex.quote(x) for x in argv)
        rc = shell.run(cmd, shell_mode=True)

        # Collect produced artifacts
        conv_fp = os.path.join(work, f"{out_prefix}_converted.csv")

        converted_tbl = None
        if os.path.exists(conv_fp) and os.path.getsize(conv_fp) > 0:
            converted_tbl = TableImporter.call(
                File(conv_fp),
                {"delimiter": ",", "header": 0, "file_format": "csv"},
            )

        # Annotated file may be CSV or TSV depending on input delimiter
        ann_csv = os.path.join(work, f"{out_prefix}_annotated.csv")
        ann_tsv = os.path.join(work, f"{out_prefix}_annotated.tsv")
        annotated_fp = ann_csv if os.path.exists(ann_csv) else (ann_tsv if os.path.exists(ann_tsv) else None)

        annotated_tbl = None
        if annotated_fp and os.path.getsize(annotated_fp) > 0:
            # pick delimiter from extension
            if annotated_fp.lower().endswith(".tsv"):
                importer_params = {"delimiter": "\t", "header": 0, "file_format": "tsv"}
            else:
                importer_params = {"delimiter": ",", "header": 0, "file_format": "csv"}
            annotated_tbl = TableImporter.call(File(annotated_fp), importer_params)

        # Return as Tables
        return {
            "converted_table": converted_tbl,
            "annotated_file": annotated_tbl,
        }