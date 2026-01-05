# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    FloatParam,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    PlotlyResource,
    ResourceSet,
    ShellProxy,
    StrParam,
    Table,
    TableExporter,
    TableImporter,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .ora_analysis_env import OraAnalysisShellProxyHelper

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
"Prunus dulcis","Pomacea canaliculata","Meloidogyne enterolobii","Halicephalobus sp. NKZ332","Halyomorpha halys","Albugo laibachii Nc14",
"Drosophila busckii","Adineta vaga","Saccharum spontaneum","Xiphophorus couchianus","Opisthorchis felineus","Meloidogyne arenaria",
"Lactuca sativa","Aspergillus flavus NRRL3357","Glossina pallidipes","Chelonus insularis","Chondrus crispus","Saccoglossus kowalevskii",
"Sorex araneus","Pyricularia oryzae 70-15","Candidozyma haemuli","Anopheles funestus","Ichneumon xanthorius","Rhagoletis pomonella",
"Drosophila pseudoobscura","Trialeurodes vaporariorum","Aphelenchoides besseyi","Amphiprion ocellaris","Glyphotaelius pellucidus",
"Ixodes persulcatus","Caenorhabditis panamensis","Salmo salar","Leishmania major strain Friedlin","Actinia equina","Panagrellus redivivus",
"Steinernema glaseri","Manihot esculenta","Leersia perrieri","Sciurus vulgaris","Nymphaea colorata","Heterodera glycines","Drosophila elegans",
"Toxocara canis","Candida parapsilosis CDC317","Solanum lycopersicum","Bactrocera dorsalis","Prunus avium","Trichinella sp. T9",
"Limnephilus marmoratus","Wuchereria bancrofti","Bison bison bison","Aspergillus fischeri NRRL 181","Micoletzkya japonica",
"Puccinia striiformis f. sp. tritici PST-130","Trypanosoma brucei","Erinaceus europaeus","Plasmodium falciparum 3D7",
"Trichinella spiralis","Bos mutus","Aspergillus clavatus NRRL 1","Clonorchis sinensis","Lineus longissimus","Anopheles christyi",
"Oryza glumipatula","Rhabditophanes sp. KR3021","Denticeps clupeoides","Galendromus occidentalis","Capsicum annuum","Phlebotomus papatasi",
"Schmidtea mediterranea","Geospiza fortis","Brugia malayi","Astyanax mexicanus","Dinothrombium tinctorium","Selaginella moellendorffii",
"Hyaloperonospora arabidopsidis Emoy2","Gigantopelta aegis","Globisporangium ultimum DAOM BR144","Pseudo-nitzschia multistriata",
"Ancylostoma ceylanicum","Myripristis murdjan","Serinus canaria","Gadus morhua","Citrullus lanatus","Opisthorchis viverrini",
"Ictidomys tridecemlineatus","Cyclopterus lumpus","Ctenocephalides felis","Electrophorus electricus","Clupea harengus",
"Populus trichocarpa","Plasmodium knowlesi strain H","Limnephilus rhombicus","Guillardia theta CCMP2712","Heliconius melpomene",
"Dibothriocephalus latus","Hucho hucho","Acrobeloides nanus","Mustela putorius furo","Nomascus leucogenys","Ustilago maydis",
"Prunus persica","Pectinophora gossypiella","Oryzias sinensis","Ooceraea biroi","Caenorhabditis remanei","Paragonimus skrjabini miyazakii",
"Meloidogyne floridensis","Harpegnathos saltator","Ascaris suum","Ovis aries","Orbicella faveolata","Globisporangium irregulare DAOM BR486",
"Bos indicus x Bos taurus","Panicum hallii","Pristionchus entomophagus","Caenorhabditis auriculariae","Phytopythium vexans DAOM BR484",
"Dendroctonus ponderosae","Gasterosteus aculeatus aculeatus","Schizosaccharomyces japonicus yFS275","Patiria miniata",
"Meloidogyne chitwoodi","Biomphalaria glabrata","Oryza rufipogon","Aedes albopictus","Dracunculus medinensis","Gyrodactylus bullatarudis",
"Medicago truncatula","Danio rerio","Carlito syrichta","Sabatieria punctata","Saimiri boliviensis boliviensis","Octopus bimaculoides",
"Panagrolaimus davidi","Microctonus hyperodae","Musa acuminata subsp. malaccensis","Leptotrombidium deliense","Bos grunniens",
"Brassica juncea","Sesamum indicum","Anopheles quadriannulatus","Oesophagostomum dentatum","Onchocerca flexuosa","Ficedula albicollis",
"Lytechinus variegatus","Equus asinus","Aspergillus fumigatus A1163","Globisporangium iwayamae DAOM BR242034","Ancylostoma duodenale",
"Phytophthora kernoviae 00238/432","Paramacrobiotus metropolitanus","Equus caballus","Microbotryum lychnidis-dioicae p1A1 Lamole",
"Zea mays","Drosophila albomicans","Phytophthora infestans T30-4","Drosophila takahashii","Phytophthora ramorum","Botrytis cinerea B05.10",
"Haemonchus contortus","Gyrodactylus salaris","Notechis scutatus","Steinernema hermaphroditum","Nannospalax galili","Crocodylus porosus",
"Dermatophagoides pteronyssinus","Brugia pahangi","Coturnix japonica","Ancylostoma caninum","Ursus maritimus","Panicum hallii var. hallii",
"Schizosaccharomyces octosporus yFS286","Penaeus vannamei","Bradynema listronoti","Drosophila navojoa","Parus major",
"Dictyostelium discoideum AX4","Macaca mulatta","Nicotiana attenuata","Arabidopsis lyrata subsp. lyrata","Heterorhabditis bacteriophora",
"Quercus suber","Exaiptasia diaphana","Stylophora pistillata","Fusarium verticillioides 7600","Dipodomys ordii","Portunus trituberculatus",
"Thrips palmi","Melampsora larici-populina 98AG31","Trichinella sp. T8","Habropoda laboriosa","Camelina sativa","Ostrea edulis",
"Pteropus vampyrus","Anopheles gambiae","Malus domestica","Bombus huntii","Tribolium madens","Drosophila grimshawi","Cottoperca gobio",
"Theristus sp. LFF4_11","Leptinotarsa decemlineata","Cimex lectularius","Camponotus floridanus","Anopheles sinensis","Rodentolepis nana",
"Melanaphis sacchari","Scophthalmus maximus","Chelonoidis abingdonii","Pyricularia oryzae","Paralinhomoeus sp. GSCO2_6",
"Callorhinchus milii","Kalanchoe fedtschenkoi","Trichoderma virens Gv29-8","Helobdella robusta","Drosophila rhopaloa",
"Puccinia graminis f. sp. tritici 04KEN156/4","Fusarium oxysporum f. sp. lycopersici 4287","Papio anubis","Myotis lucifugus",
"Setaria italica","Dicrocoelium dendriticum","Olea europaea subsp. europaea","Gaeumannomyces tritici R3-111a-1",
"Fusarium pseudograminearum CS3096","Tupaia belangeri","Sitophilus oryzae","Pythium arrhenomanes ATCC 12531","Pocillopora damicornis",
"Haemaphysalis longicornis","Apis dorsata","Haliotis rubra","Microlaimidae sp. YZB2_3","Strigops habroptila",
"Stegodyphus dumicola","Protopolystoma xenopodis","Anopheles culicifacies","Dasypus novemcinctus","Maylandia zebra",
"Mesodorylaimus sp. YZB2_4","Strongyloides venezuelensis","Brugia timori","Gopherus evgoodei","Candidozyma pseudohaemuli",
"Brassica rapa","Triticum dicoccoides","Glycine soja","Acanthochromis polyacanthus","Pogonomyrmex barbatus","Parascaris univalens",
"Pristionchus pacificus","Bursaphelenchus okinawaensis","Mastacembelus armatus","Arabidopsis halleri subsp. gemmifera",
"Podarcis muralis","Amblyteles armatorius","Plectus sambesii","Gallus gallus","Mus caroli","Eucalyptus grandis","Procambarus clarkii",
"Drosophila gunungcola","Notamacropus eugenii","Anser brachyrhynchus","Zootermopsis nevadensis","Prolemur simus",
"Echinococcus multilocularis","Amphiprion percula","Octopus sinensis","Drosophila willistoni","Cercocebus atys","Clytia hemisphaerica",
"Dirofilaria immitis","Aspergillus oryzae RIB40","Rhinopithecus roxellana","Aphanomyces invadans","Poecilia formosa","Drosophila kikkawai",
"Homo sapiens","Homalodisca vitripennis","Fusarium culmorum CS7071","Aphanomyces astaci","Panagrolaimus superbus",
"Pelodiscus sinensis","Blumeria hordei DH14","Betta splendens","Daphnia pulicaria","Sander lucioperca","Glossina fuscipes","Anopheles farauti",
"Asparagus officinalis","Cricetulus griseus","Acanthocheilonema viteae","Candidozyma auris","Penaeus japonicus",
"Drosophila miranda","Brassica oleracea var. oleracea","Propithecus coquereli","Hyalomma asiaticum","Urocitellus parryii",
"Aspergillus fumigatus Af293","Giardia lamblia ATCC 50803","Anopheles maculatus","Triticum turgidum subsp. durum","Chara braunii",
"Cryptotermes secundus","Parascaris equorum","Coffea canephora","Aphelenchoides bicaudatus","Ficus carica","Vulpes vulpes",
"Anopheles stephensi","Anisakis simplex","Rhipicephalus microplus","Oryza sativa aromatic subgroup","Histoplasma ohiense (nom. inval.)",
"Hypsibius exemplaris","Verticillium dahliae JR2","Puccinia triticina 1-1 BBBD Race 1","Phascolarctos cinereus","Pecten maximus",
"Oppia nitens","Esox lucius","Bactrocera neohumeralis","Panagrolaimus sp. PS1159","Dermacentor silvarum","Nippostrongylus brasiliensis",
"Magnaporthiopsis poae ATCC 64411","Leptobrachium leishanense","Moschus moschiferus","Pomphorhynchus laevis",
"Uloborus diversus","Cervus hanglu yarkandensis","Ursus americanus","Drosophila pseudoobscura pseudoobscura","Aedes aegypti",
"Allodiplogaster sudhausi","Amphibalanus amphitrite","Brassica rapa subsp. trilocularis","Rosa chinensis","Aegilops umbellulata",
"Emiliania huxleyi CCMP1516","Trichoplax adhaerens","Anopheles coluzzii","Olea europaea var. sylvestris","Ciona savignyi",
"Polistes dominula","Stegastes partitus","Mnemiopsis leidyi","Oryctolagus cuniculus","Taenia saginata","Komagataella phaffii GS115",
"Corchorus capsularis","Tribolium castaneum","Anopheles merus","Sclerotinia sclerotiorum 1980 UF-70","Stomoxys calcitrans",
"Kryptolebias marmoratus","Caenorhabditis sinica","Callithrix jacchus","Plasmodium berghei ANKA","Daphnia carinata",
"Colletotrichum fructicola Nara gc5","Trichinella sp. T6","Caenorhabditis brenneri","Neodiprion lecontei","Panthera leo",
"Caenorhabditis latens","Cataglyphis hispanica","Choloepus hoffmanni","Labrus bergylta","Plenodomus lingam JN3",
"Corymbia citriodora subsp. variegata","Epsilonema sp. ZAB3_2","Nilaparvata lugens","Sipha flava","Papaver somniferum",
"Sparus aurata","Oryza brachyantha","Procavia capensis","Pristionchus fissidentatus","Drosophila melanogaster","Monodelphis domestica",
"Ipomoea triloba","Patella pellucida","Pygocentrus nattereri","Meloidogyne hapla","Triticum urartu","Ornithorhynchus anatinus",
"Toxoplasma gondii ME49","Microcebus murinus","Meloidogyne javanica","Aplysia californica","Glycine max","Chinchilla lanigera",
"Phaeodactylum tricornutum CCAP 1055/1","Penaeus monodon","Drosophila erecta","Sarcophilus harrisii","Caenorhabditis angaria",
"Gongylonema pulchrum","Polistes canadensis","Pristionchus japonicus","Schistocerca americana","Hymenolepis diminuta",
"Ochotona princeps","Angiostrongylus costaricensis","Schizosaccharomyces cryophilus OY26","Haemonchus placei","Atta cephalotes",
"Rhinopithecus bieti","Delphinapterus leucas","Monomorium pharaonis","Amyelois transitella","Candida albicans SC5314","Oryza barthii",
"Belgica antarctica","Mus musculus","Taeniopygia guttata","Trichoderma reesei QM6a","Magallana gigas","Balaenoptera musculus",
"Aotus nancymaae","Mesorhabditis belari","Diaphorina citri","Lates calcarifer","Nothobranchius furzeri","Drosophila subobscura",
"Daphnia pulex","Hymenolepis microstoma","Pistacia vera","Oryza punctata","Cyanidioschyzon merolae strain 10D","Drosophila suzukii",
"Mandrillus leucophaeus","Drosophila ananassae","Trichinella britovi","Bombus impatiens","Colletotrichum graminicola M1.001",
"Echinochloa crus-galli","Steinernema monticolum","Anopheles albimanus","Anabas testudineus","Trileptium ribeirensis",
"Colletotrichum orbiculare MAFF 240422","Trifolium pratense","Capra hircus","Capitella teleta","Oscheius tipulae","Diploscapter coronatus",
"Loxodonta africana","Crassostrea virginica","Panthera tigris altaica","Trichuris trichiura","Xiphophorus maculatus",
"Hofstenia miamia","Octodon degus","Schistosoma mansoni","Schistocerca cancellata","Zerene cesonia","Megaselia scalaris",
"Venturia canescens","Meleagris gallopavo","Eragrostis curvula","Camelus dromedarius","Drosophila subpulchrella","Bactrocera tryoni",
"Bunonema sp. RGD898","Schistosoma guineensis","Monodon monoceros","Actinidia chinensis var. chinensis","Megachile rotundata",
"Rhipicephalus sanguineus","Caenorhabditis sulstoni","Candida tropicalis MYA-3404","Citrus x clementina","Orchesella cincta",
"Mus spicilegus","Limulus polyphemus","Limnoperna fortunei","Leguminivora glycinivorella","Sphenodon punctatus","Folsomia candida",
"Drosophila santomea","Cucumis sativus","Colletotrichum higginsianum","Strigamia maritima","Scleropages formosus","Neodiprion pinetum",
"Culex quinquefasciatus","Chrysemys picta bellii","Asterias rubens","Oryza meridionalis","Ceratitis capitata","Pongo abelii",
"Neurospora crassa OR74A","Orussus abietinus","Drosophila mauritiana","Ptycholaimellus sp. GST1_10","Schistosoma spindalis",
"Strongyloides stercoralis","Drosophila teissieri","Trichobilharzia szidati","Enterobius vermicularis","Neolamprologus brichardi",
"Setaria digitata","Bombus vancouverensis nearcticus","Trichinella zimbabwensis","Cucumis melo","Apis florea","Hydractinia symbiolongicarpus",
"Thelohanellus kitauei","Schistocephalus solidus","Pseudonaja textilis","Pan paniscus","Chlorocebus sabaeus","Fundulus heteroclitus",
"Phytophthora nicotianae P1569","Branchiostoma lanceolatum","Mesocestoides corti","Pyrenophora tritici-repentis Pt-1C-BFP",
"Peromyscus maniculatus bairdii","Dufourea novaeangliae","Mesorhabditis spiculigera","Yarrowia lipolytica CLIB122","Fasciola gigantica",
"Vigna angularis","Trichinella papuae","Soboliphyme baturini","Mya arenaria","Cotesia glomerata","Aegilops tauschii subsp. strangulata",
"Linepithema humile","Seriola lalandi dorsalis","Homarus americanus","Schizosaccharomyces pombe 972h-","Amphilophus citrinellus",
"Oryzias javanicus","Dermacentor andersoni","Pisum sativum","Astatotilapia calliptera","Glossina austeni","Gorilla gorilla gorilla",
"Hydatigera taeniaeformis","Anopheles minimus","Fraxinus excelsior","Caenorhabditis japonica","Cylicostephanus goldi","Polistes fuscatus",
"Steinernema scapterisci","Anopheles darlingi","Helicoverpa armigera","Caenorhabditis becei","Penaeus chinensis","Pundamilia nyererei",
"Drosophila obscura","Puccinia graminis f. sp. tritici CRL 75-36-700-3","Latimeria chalumnae","Anopheles epiroticus","Phlebotomus perniciosus",
"Parasteatoda tepidariorum","Culicoides sonorensis","Fusarium vanettenii 77-13-4","Chlamydomonas reinhardtii","Homarus gammarus",
"Echinococcus canadensis","Taenia solium","Chenopodium quinoa","Oryza nivara","Patella vulgata","Parelaphostrongylus tenuis",
"Koerneria luziae","Mus pahari","Tuber melanosporum Mel28","Daphnia magna","Bos taurus","Pristionchus exspectatus",
"Bigelowiella natans CCMP2755","Eufriesea mexicana","Ostreococcus lucimarinus CCE9901","Copidosoma floridanum","Tigriopus californicus",
"Echinococcus oligarthrus","Schistocerca piceifrons","Eurytemora affinis","Triticum spelta","Echinops telfairi","Trissonchulus latispiculum",
"Aspergillus niger","Pythium aphanidermatum DAOM BR444","Onthophagus taurus","Ancistrocerus nigricornis","Bombyx mori","Oncorhynchus kisutch",
"Fasciolopsis buskii","Eptatretus burgeri","Drosophila bipectinata","Acromyrmex echinatior","Onchocerca volvulus","Vigna unguiculata",
"Rattus norvegicus","Romanomermis culicivorax","Caenorhabditis bovis","Sporisorium reilianum SRZ2"
]

def _grid_from_long(
    long_df: pd.DataFrame,
    title: str,
    genes_per_page: int = 35,
    terms_per_page: int = 35
) -> PlotlyResource:
    """
    Interactive Terms×Genes grid with:
      - Full Y labels: "term_name | term_id | source | p-value"
      - Vertical paging for terms (dropdown), generous line spacing
      - Genes pagination (dropdown)
      - Extra y-range padding to avoid clipping of bubbles
      - Wider canvas and larger margins for breathing room
    Always returns a valid PlotlyResource (empty figure if long_df is empty).
    """
    if long_df is None or long_df.empty:
        fig = go.Figure()
        fig.update_layout(
            title=title,
            width=1600, height=480,
            margin=dict(l=140, r=90, t=130, b=170),
            xaxis_title="<b>Genes</b>",
            yaxis_title="<b>Terms</b>",
        )
        return PlotlyResource(fig)

    d = long_df.copy()
    src_order = {"GO:CC": 0, "GO:BP": 1, "GO:MF": 2, "KEGG": 3}
    d["_src_ord"] = d["source"].map(src_order).fillna(9)
    d["_term_key"] = d["term_id"].astype(str) + "||" + d["source"].astype(str)

    # Order terms by source then score
    best_score = d.groupby("_term_key")["score"].max()
    first_rows = d.drop_duplicates("_term_key").set_index("_term_key")

    def fmt_p(x):
        try:
            return f"{float(x):.2e}"
        except Exception:
            return "NA"

    y_labels_full = (
        first_rows["term_name"].astype(str) + " | " +
        first_rows["term_id"].astype(str)   + " | " +
        first_rows["source"].astype(str)    + " | " +
        first_rows["pvalue"].map(fmt_p)
    )

    order = (
        pd.DataFrame({
            "key": first_rows.index,
            "src": first_rows["_src_ord"].values,
            "score": best_score[first_rows.index].values,
            "label": y_labels_full.values,
        })
        .sort_values(["src", "score"], ascending=[True, False])
        .reset_index(drop=False)
    )
    n_terms = len(order)

    # Interline generous
    per_row_px = int(np.interp(n_terms, [0, 30, 60, 90, 120], [56, 48, 40, 34, 30]))
    total_height = int(220 + per_row_px * max(1, n_terms))
    total_height = min(total_height, 8000)

    # Y positions & label lookup
    y_positions = np.arange(n_terms, dtype=float)
    key_to_y: Dict[str, float] = dict(zip(order["key"], y_positions))
    term_lookup = dict(zip(order["key"], order["label"]))

    # Color by p-value (reversed Viridis)
    pvals = pd.to_numeric(d["pvalue"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    max95 = np.nanpercentile(pvals, 95) if np.isfinite(pvals).any() else 1.0
    d["_p"] = np.clip(pvals, 0, max95)

    # Bubble size by |log2FC|
    abs_lfc = pd.to_numeric(d.get("log2fc", np.nan), errors="coerce").abs()
    lfc95   = np.nanpercentile(abs_lfc, 95) if np.isfinite(abs_lfc).any() else 1.0
    scale   = (abs_lfc / (lfc95 if lfc95 > 0 else 1.0)).fillna(0.0)
    min_bub = max(7.0, per_row_px * 0.28)
    max_bub = max(min_bub + 5.0, per_row_px * 0.58)
    d["_sz"] = (min_bub + (max_bub - min_bub) * np.clip(scale, 0, 1)).astype(float)

    # Term×gene presence
    presence = d[["_term_key","gene_label","_p","_sz","log2fc"]].drop_duplicates()
    presence["_y"] = presence["_term_key"].map(key_to_y)

    # Pagination (X) by genes
    genes = sorted(presence["gene_label"].dropna().astype(str).unique().tolist())
    pages_g: List[List[str]] = [genes[i:i+genes_per_page] for i in range(0, len(genes), genes_per_page)] or [genes]

    # One trace per gene page; vertical paging via yaxis windows
    fig = go.Figure()
    for gi, page_genes in enumerate(pages_g, start=1):
        sub = presence[presence["gene_label"].isin(page_genes)]
        xmap = {g: i for i, g in enumerate(page_genes)}
        x = [xmap[g] for g in sub["gene_label"]]
        y = sub["_y"].astype(float).values

        fig.add_trace(go.Scatter(
            x=x, y=y, mode="markers",
            marker=dict(
                size=sub["_sz"], color=sub["_p"],
                colorscale="Viridis", reversescale=True,
                colorbar=dict(title="p-value", thickness=12, len=0.85),
                line=dict(width=0.3, color="#444"),
                opacity=0.9,
            ),
            cliponaxis=False,  # avoid clipping at edges
            visible=(gi == 1),
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "Gene: %{customdata[1]}<br>"
                "log2FC: %{customdata[2]}<br>"
                "p-value: %{customdata[3]}<extra></extra>"
            ),
            customdata=np.stack([
                sub["_term_key"].map(term_lookup).values,
                sub["gene_label"].astype(str).values,
                pd.to_numeric(sub["log2fc"], errors="coerce").round(3).astype(object).values,
                sub["_p"].astype(float).map(lambda v: f"{v:.2e}").values
            ], axis=1),
            name=f"Genes page {gi}",
        ))

    # Dropdown 1: genes pagination (switch trace + update x ticks)
    buttons_genes = []
    for gi, page_genes in enumerate(pages_g, start=1):
        vis = [False]*len(pages_g); vis[gi-1] = True
        buttons_genes.append(dict(
            label=f"Genes {gi}/{len(pages_g)}",
            method="update",
            args=[{"visible": vis},
                  {"xaxis": {
                      "tickmode": "array",
                      "tickvals": list(range(len(page_genes))),
                      "ticktext": page_genes,
                      "tickangle": 45,
                      "tickfont": {"size": 10}
                  }}]
        ))

    # Dropdown 2: vertical paging by terms
    tp = max(1, int(terms_per_page))
    term_windows = [(s, min(s+tp, n_terms)) for s in range(0, n_terms, tp)]
    ypad = 0.7  # padding in "row units" to avoid top/bottom clipping

    def _y_layout_for_window(win):
        s, e = win
        tickvals = y_positions[s:e].tolist()
        ticktext = order.loc[s:e-1, "label"].tolist()
        return dict(
            tickmode="array",
            tickvals=tickvals,
            ticktext=ticktext,
            tickfont=dict(size=11 if tp <= 40 else 10),
            autorange=False,
            range=[(e - 0.5) + ypad, (s - 0.5) - ypad],  # reversed axis with padding
            domain=[0.12, 0.98],  # extra space above/below plotting area
            showline=False, zeroline=False,
            showgrid=True, gridcolor="#f7f7f7", gridwidth=0.6,
        )

    buttons_terms = []
    for ti, win in enumerate(term_windows, start=1):
        buttons_terms.append(dict(
            label=f"Terms {ti}/{len(term_windows)}",
            method="relayout",
            args=[{"yaxis": _y_layout_for_window(win)}]
        ))

    # Adaptive left margin to fit long labels
    max_lab_len = int(order["label"].str.len().max())
    left_margin = min(140 + 7 * max_lab_len, 1600)

    # Initial layout
    fig.update_layout(
        title=title,
        font=dict(size=10.5),
        title_font=dict(size=20),
        xaxis=dict(
            title="<b>Genes</b>",
            title_standoff=16,
            title_font=dict(size=13, family="Arial"),
            tickmode="array",
            tickvals=list(range(len(pages_g[0]))),
            ticktext=pages_g[0],
            tickangle=45,
            tickfont=dict(size=10),
            showline=False, zeroline=False,
            showgrid=True, gridcolor="#f0f0f0", gridwidth=0.6,
        ),
        yaxis=_y_layout_for_window(term_windows[0]),
        updatemenus=[
            dict(type="dropdown", x=0.80, y=1.13, showactive=True,
                 direction="down", buttons=buttons_terms),
            dict(type="dropdown", x=0.90, y=1.13, showactive=True,
                 direction="down", buttons=buttons_genes),
        ],
        margin=dict(l=left_margin, r=90, t=140, b=180),
        width=1650,
        height=min(total_height + 120, 8200),
    )

    # Header above labels with extra gap to the first row
    fig.add_annotation(
        xref="paper", yref="paper",
        x=-0.40, y=1.07,
        xanchor="left",
        text="<b>term_name | term_id | source | p-value</b>",
        showarrow=False, font=dict(size=12)
    )

    # If very few genes, lightly tighten X range
    if len(pages_g[0]) <= 5:
        fig.update_xaxes(range=[-0.6, max(1.6, len(pages_g[0]) - 0.4)])

    return PlotlyResource(fig)

# ──────────────────────────── Task runner ───────────────────────────────
@task_decorator(
    "ORAEnrichmentMultiSpecies",
    human_name="Functional enrichment analysis based on ORA",
    short_description="g:Profiler ORA of DE genes — tables, PNG barplots, and an interactive Terms×Genes grid.",
)
class ORAEnrichmentTask(Task):
    """
    This task perform functional enrichment (ORA) of differentially expressed genes using g:Profiler,
    then package the results as tables, static barplots, and an interactive Plotly “Terms × Genes” grid.
    The analysis script reads a DE results CSV or Table, builds three gene sets (ALL, UP, DOWN) based on FDR (padj) and optional |log2FC| thresholds, runs g:Profiler for selected sources (GO BP/MF/CC, KEGG), filters by FDR, writes CSVs, and saves top-N barplots.
    The runner script executes that analysis, imports all outputs, and builds an interactive grid where bubble color reflects
    p-value and size reflects |log2FC|, with dropdown pagination for both genes (horizontal) and terms (vertical)
    via configurable grid_genes_per_page and grid_terms_per_page, plus extra spacing and margins to keep labels readable.
    """

    input_specs: Final[InputSpecs] = InputSpecs({
        "de_table_file": InputSpec((Table, File), human_name="DE/Full results (Table or CSV file)"),
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "tables":         OutputSpec(ResourceSet, human_name="ORA tables"),
        "barplots":       OutputSpec(ResourceSet, human_name="Barplots (top-N)"),
        "term_gene_grid": OutputSpec(ResourceSet, human_name="Interactive Terms×Genes"),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "organism_name": StrParam(allowed_values=SCIENTIFIC_NAMES,
            short_description="Scientific name (select from list)."),
        "genes_colname": StrParam(
            short_description="Name of the column that holds gene IDs (Ensembl or symbols or mixed)."),
        "padj_threshold": FloatParam(default_value=0.05, min_value=0.0,
            short_description="FDR threshold used to filter g:Profiler results."),
        "abs_log2fc": FloatParam(default_value=0.0,
            short_description="Absolute |log2FC| (used only upstream; not by g:Profiler)."),
        "topn_plot": FloatParam(default_value=20.0, min_value=5.0,
            short_description="Top-N terms for static barplots."),
        "sources_list": StrParam(default_value="GO:BP,GO:MF,GO:CC,KEGG",
            short_description="Comma-separated sources for g:Profiler."),
        "grid_genes_per_page": FloatParam(default_value=25.0, min_value=1.0,
            short_description="Genes per page in the interactive grid."),
        "grid_terms_per_page": FloatParam(default_value=20.0, min_value=1.0,
            short_description="Terms per page in the interactive grid."),
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_ora_analysis.py",  # the analysis script
    )

    def run(self, p: ConfigParams, ins: TaskInputs) -> TaskOutputs:
        # ── Accept Table or File for DE results; export Table -> CSV if needed ──
        de_res = ins["de_table_file"]
        if isinstance(de_res, Table):
            de_file: File = TableExporter.call(de_res, {"file_format": "csv"})
            de_csv = de_file.path
        elif isinstance(de_res, File):
            de_csv = de_res.path
        else:
            raise RuntimeError("Unsupported type for 'de_table_file' (expected Table or File).")

        species = p["organism_name"].strip()
        idcol   = p["genes_colname"].strip()
        padj    = float(p["padj_threshold"])
        lfc     = float(p["abs_log2fc"])
        topn    = int(p["topn_plot"])
        sources = p["sources_list"].strip() or "GO:BP,GO:MF,GO:CC,KEGG"
        genes_per_page = int(p.get("grid_genes_per_page", 35))
        terms_per_page = int(p.get("grid_terms_per_page", 35))

        shell: ShellProxy = OraAnalysisShellProxyHelper.create_proxy(self.message_dispatcher)
        work = shell.working_dir
        csv_out = os.path.join(work, "csv")
        fig_out = os.path.join(work, "figs")

        cmd = (
            'python3 {py} --infile "{de}" --species "{sp}" '
            '--id_column "{idc}" --padj_thr {padj} --lfc_thr {lfc} '
            '--sources "{src}" --topn {topn} '
            '--outdir "{outdir}" --csv_dir "{csvdir}" --fig_dir "{figdir}"'
        ).format(
            py=self.python_file_path, de=de_csv, sp=species.replace('"', '\\"'),
            idc=idcol, padj=padj, lfc=lfc, src=sources, topn=topn,
            outdir=work, csvdir=csv_out, figdir=fig_out,
        )

        if shell.run(cmd, shell_mode=True) != 0:
            return {"tables": ResourceSet(), "barplots": ResourceSet(), "term_gene_grid": ResourceSet()}

        tables, barplots, grids = ResourceSet(), ResourceSet(), ResourceSet()

        # Collect tables and remember all long-matrix CSVs (even empty) for grids
        long_csvs: List[str] = []
        if os.path.isdir(csv_out):
            for fn in sorted(os.listdir(csv_out)):
                fp = os.path.join(csv_out, fn)
                if fn.endswith("_matrix_long.csv"):
                    long_csvs.append(fp)  # include empty; we still create a grid
                elif fn.lower().endswith(".csv"):
                    try:
                        if os.path.getsize(fp) > 0:
                            tables.add_resource(
                                TableImporter.call(
                                    File(fp),
                                    {"delimiter": ",", "header": 0, "file_format": "csv"},
                                ),
                                fn,
                            )
                    except Exception:
                        pass

        # Collect barplots
        if os.path.isdir(fig_out):
            for fn in sorted(os.listdir(fig_out)):
                if fn.lower().endswith(".png"):
                    fp = os.path.join(fig_out, fn)
                    try:
                        barplots.add_resource(File(fp), Path(fn).stem)
                    except Exception:
                        pass

        # Build grids for all *_matrix_long.csv
        for fp in long_csvs:
            try:
                d = pd.read_csv(fp) if os.path.getsize(fp) > 0 else pd.DataFrame()
            except Exception:
                d = pd.DataFrame()
            title = f"Terms × Genes — {Path(fp).stem.replace('_matrix_long','')}"
            try:
                pr = _grid_from_long(
                    d, title=title,
                    genes_per_page=genes_per_page,
                    terms_per_page=terms_per_page
                )
                grids.add_resource(pr, Path(fp).stem)
            except Exception:
                # Fallback: empty figure with title
                fig = go.Figure()
                fig.update_layout(title=title, width=1200, height=400)
                grids.add_resource(PlotlyResource(fig), Path(fp).stem)

        return {"tables": tables, "barplots": barplots, "term_gene_grid": grids}
