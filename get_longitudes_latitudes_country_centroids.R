library(rgeos)
library(rworldmap)

# get world map
wmap <- getMap(resolution="high")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with longitudes & latitudes of country centroids
countries_long_lat <- as.data.frame(centroids)
colnames(countries_long_lat) = c("longitude","latitude")
countries_long_lat$country = rownames(countries_long_lat)
rownames(countries_long_lat) = NULL
write.csv(countries_long_lat, ".//data//GISAID//countries_long_lat_centroids.csv", row.names=F)
countries_long_lat
#        longitude    latitude                                  country
# 1    -69.9826747  12.5208888                                    Aruba
# 2     66.0047311  33.8352322                              Afghanistan
# 3     17.5373615 -12.2933611                                   Angola
# 4    -63.0649846  18.2239669                                 Anguilla
# 5     20.0498303  41.1424513                                  Albania
# 6     19.9532464  60.2149004                                    Aland
# 7      1.5605337  42.5422917                                  Andorra
# 8     54.3001574  23.9052718                     United Arab Emirates
# 9    -65.1798088 -35.3813535                                Argentina
# 10    44.9299343  40.2895221                                  Armenia
# 11  -170.7180186 -14.3044892                           American Samoa
# 12    19.9211327 -80.5085641                               Antarctica
# 13   123.5838304 -12.4299425              Ashmore and Cartier Islands
# 14    69.2267342 -49.2489599      French Southern and Antarctic Lands
# 15   -61.7946956  17.2775426                      Antigua and Barbuda
# 16   134.4909988 -25.7328874                                Australia
# 17    14.1264781  47.5854985                                  Austria
# 18    47.5459965  40.2882767                               Azerbaijan
# 19    29.8751201  -3.3593926                                  Burundi
# 20     4.6406457  50.6398144                                  Belgium
# 21     2.3278453   9.6417533                                    Benin
# 22    -1.7545684  12.2695368                             Burkina Faso
# 23    90.2381279  23.8673041                               Bangladesh
# 24    25.2155233  42.7688999                                 Bulgaria
# 25    50.5419776  26.0420085                                  Bahrain
# 26   -76.6283574  24.2903122                              The Bahamas
# 27    17.7687659  44.1744991                   Bosnia and Herzegovina
# 28   -62.8406871  17.8988104                         Saint Barthelemy
# 29    28.0320871  53.5313115                                  Belarus
# 30   -88.7101023  17.2002747                                   Belize
# 31   -64.7545418  32.3136859                                  Bermuda
# 32   -64.6853892 -16.7081454                                  Bolivia
# 33   -53.0978320 -10.7877757                                   Brazil
# 34   -59.5598020  13.1814537                                 Barbados
# 35   114.7220255   4.5196824                                   Brunei
# 36    90.4018950  27.4110619                                   Bhutan
# 37    23.7985373 -22.1840203                                 Botswana
# 38    20.4682636   6.5682300                 Central African Republic
# 39   -98.3077572  61.3620683                                   Canada
# 40     8.2086601  46.7978530                              Switzerland
# 41   -71.3825588 -37.7306975                                    Chile
# 42   103.8190738  36.5617680                                    China
# 43    -5.5692160   7.6284287                              Ivory Coast
# 44    12.7396427   5.6911030                                 Cameroon
# 45    23.6439564  -2.8774601         Democratic Republic of the Congo
# 46    15.2196572  -0.8378817                    Republic of the Congo
# 47  -159.7872415 -21.2192681                             Cook Islands
# 48   -73.0811494   3.9138316                                 Colombia
# 49    43.6823810 -11.8777757                                  Comoros
# 50   -23.9599097  15.9551882                               Cape Verde
# 51   -84.1920879   9.9763424                               Costa Rica
# 52   -79.0160691  21.6229002                                     Cuba
# 53   -68.9712225  12.1955244                                  Curacao
# 54   -80.9121306  19.4289645                           Cayman Islands
# 55    33.5684401  35.2627605                          Northern Cyprus
# 56    33.0060036  34.9166693                                   Cyprus
# 57    15.3124007  49.7334107                           Czech Republic
# 58    10.3857769  51.1069790                                  Germany
# 59    42.5606768  11.7487205                                 Djibouti
# 60   -61.3577312  15.4394721                                 Dominica
# 61    10.0279942  55.9812639                                  Denmark
# 62   -70.5056909  18.8943308                       Dominican Republic
# 63     2.6173210  28.1589384                                  Algeria
# 64   -78.7520357  -1.4238195                                  Ecuador
# 65    29.8618993  26.4959311                                    Egypt
# 66    38.8461544  15.3618688                                  Eritrea
# 67    -3.6475487  40.2444863                                    Spain
# 68    25.5424855  58.6719242                                  Estonia
# 69    39.6008014   8.6227842                                 Ethiopia
# 70    26.2746702  64.4988484                                  Finland
# 71   165.4505548 -17.4285609                                     Fiji
# 72   -59.3524005 -51.7448299                         Falkland Islands
# 73     2.5361847  46.1870058                                   France
# 74    -6.8809632  62.0538623                            Faroe Islands
# 75   153.2370611   7.4528106           Federated States of Micronesia
# 76    11.7886301  -0.5866033                                    Gabon
# 77    34.3475712  31.3915035                                     Gaza
# 78    -2.8656336  54.1238716                           United Kingdom
# 79    43.5077900  42.1685594                                  Georgia
# 80    -2.5724016  49.4680986                                 Guernsey
# 81    -1.2167690   7.9534512                                    Ghana
# 82   -10.9406638  10.4362159                                   Guinea
# 83   -15.3960181  13.4496538                                   Gambia
# 84   -14.9497294  12.0474445                            Guinea Bissau
# 85    10.3413890   1.7055404                        Equatorial Guinea
# 86    22.9555621  39.0746697                                   Greece
# 87   -61.6821942  12.1172555                                  Grenada
# 88   -41.3419103  74.7105137                                Greenland
# 89   -90.3648229  15.6940416                                Guatemala
# 90   144.7679224  13.4416655                                     Guam
# 91   -58.9820260   4.7937805                                   Guyana
# 92   114.1137986  22.3982878                         Hong Kong S.A.R.
# 93    73.5205271 -53.0872530        Heard Island and McDonald Islands
# 94   -86.6151685  14.8268823                                 Honduras
# 95    16.4041308  45.0804728                                  Croatia
# 96   -72.6852736  18.9350096                                    Haiti
# 97    19.3955938  47.1627771                                  Hungary
# 98   117.2400976  -2.2150502                                Indonesia
# 99    -4.5387382  54.2241965                              Isle of Man
# 100   79.6119736  22.8857804                                    India
# 101  104.8494652 -10.6482963                 Indian Ocean Territories
# 102   72.4454024  -7.3305825           British Indian Ocean Territory
# 103   -8.1379317  53.1754426                                  Ireland
# 104   54.2740707  32.5750329                                     Iran
# 105   43.7435320  33.0397051                                     Iraq
# 106  -18.5739709  64.9957526                                  Iceland
# 107   35.0044500  31.4611019                                   Israel
# 108   12.0700095  42.7966357                                    Italy
# 109  -77.3148256  18.1569548                                  Jamaica
# 110   -2.1268765  49.2183749                                   Jersey
# 111   36.7713725  31.2457939                                   Jordan
# 112  138.0308921  37.5923006                                    Japan
# 113   77.1801129  35.3923662                          Siachen Glacier
# 114   67.2914966  48.1568795                               Kazakhstan
# 115   37.7959444   0.5998746                                    Kenya
# 116   74.5416555  41.4622179                               Kyrgyzstan
# 117  104.9069431  12.7200493                                 Cambodia
# 118  -45.6220913   0.8605225                                 Kiribati
# 119  -62.6875424  17.2645859                    Saint Kitts and Nevis
# 120  127.8391693  36.3852463                              South Korea
# 121   20.8724961  42.5707842                                   Kosovo
# 122   47.5869936  29.3343168                                   Kuwait
# 123  103.7377213  18.5021797                                     Laos
# 124   35.8801633  33.9230699                                  Lebanon
# 125   -9.3220778   6.4527827                                  Liberia
# 126   18.0086578  27.0309425                                    Libya
# 127  -60.9696933  13.8947837                              Saint Lucia
# 128    9.5357359  47.1366545                            Liechtenstein
# 129   80.7010779   7.6126631                                Sri Lanka
# 130   28.2272277 -29.5800467                                  Lesotho
# 131   23.8872004  55.3261108                                Lithuania
# 132    6.0718268  49.7672469                               Luxembourg
# 133   24.9123704  56.8508531                                   Latvia
# 134  113.5093310  22.2231081                              Macau S.A.R
# 135  -63.0597162  18.0888933                             Saint Martin
# 136   -8.4561605  29.8376279                                  Morocco
# 137    7.4062805  43.7527611                                   Monaco
# 138   28.4567374  47.1949932                                  Moldova
# 139   46.7047375 -19.3718960                               Madagascar
# 140   73.4571108   3.7285473                                 Maldives
# 141 -102.5234504  23.9475393                                   Mexico
# 142  170.3396364   7.0035280                         Marshall Islands
# 143   21.6821134  41.5953106                                Macedonia
# 144   -3.5426879  17.3458146                                     Mali
# 145   14.4052204  35.9215132                                    Malta
# 146   96.4884321  21.1856651                                  Myanmar
# 147   19.2388382  42.7889115                               Montenegro
# 148  103.0529965  46.8268143                                 Mongolia
# 149  145.6196188  15.8287301                 Northern Mariana Islands
# 150   35.5336735 -17.2738190                               Mozambique
# 151  -10.3477974  20.2573693                               Mauritania
# 152  -62.1851981  16.7394136                               Montserrat
# 153   57.5712076 -20.2777156                                Mauritius
# 154   34.2893585 -13.2180697                                   Malawi
# 155  109.6976155   3.7898726                                 Malaysia
# 156   17.2096277 -22.1303079                                  Namibia
# 157  165.6849208 -21.2999021                            New Caledonia
# 158    9.3854645  17.4191260                                    Niger
# 159  167.9492020 -29.0514384                           Norfolk Island
# 160    8.0894424   9.5941181                                  Nigeria
# 161  -85.0305289  12.8470963                                Nicaragua
# 162 -169.8699537 -19.0494623                                     Niue
# 163    5.2814735  52.1008080                              Netherlands
# 164   15.3483297  68.7501441                                   Norway
# 165   83.9158367  28.2489136                                    Nepal
# 166  166.9325852  -0.5191420                                    Nauru
# 167  171.4849539 -41.8111283                              New Zealand
# 168   56.0916597  20.6051541                                     Oman
# 169   69.3395803  29.9497522                                 Pakistan
# 170  -80.1191515   8.5175119                                   Panama
# 171 -128.3170409 -24.3649860                         Pitcairn Islands
# 172  -74.3824289  -9.1528028                                     Peru
# 173  122.8839215  11.7753974                              Philippines
# 174  134.4082215   7.2876099                                    Palau
# 175  145.2074492  -6.4641628                         Papua New Guinea
# 176   19.3901184  52.1275954                                   Poland
# 177  -66.4730731  18.2281271                              Puerto Rico
# 178  127.1924781  40.1535029                              North Korea
# 179   -8.5010551  39.5955025                                 Portugal
# 180  -58.4001408 -23.2282390                                 Paraguay
# 181 -144.9048876 -14.7222807                         French Polynesia
# 182   51.1847851  25.3060094                                    Qatar
# 183   24.9729329  45.8524294                                  Romania
# 184   96.6865735  61.9805217                                   Russia
# 185   29.9198911  -1.9903320                                   Rwanda
# 186  -12.2198297  24.2295636                           Western Sahara
# 187   44.5368704  24.1224547                             Saudi Arabia
# 188   29.9404723  15.9903543                                    Sudan
# 189   30.2479006   7.3087758                              South Sudan
# 190  -14.4734890  14.3662438                                  Senegal
# 191  103.8172532   1.3587707                                Singapore
# 192  -36.4332045 -54.4648715 South Georgia and South Sandwich Islands
# 193   -9.5476868 -12.4036762                             Saint Helena
# 194  159.6328471  -8.9217570                          Solomon Islands
# 195  -11.7927121   8.5632907                             Sierra Leone
# 196  -88.8716446  13.7394401                              El Salvador
# 197   12.4592303  43.9418694                               San Marino
# 198   46.2519833   9.7334580                               Somaliland
# 199   45.7071624   4.7506519                                  Somalia
# 200  -56.3031934  46.9191999                Saint Pierre and Miquelon
# 201   20.7895812  44.2215031                       Republic of Serbia
# 202    6.7242964   0.4439085                    Sao Tome and Principe
# 203  -55.9123483   4.1305538                                 Suriname
# 204   19.4790521  48.7054718                                 Slovakia
# 205   14.8044464  46.1155544                                 Slovenia
# 206   16.7437122  62.7923202                                   Sweden
# 207   31.4819316 -26.5584273                                Swaziland
# 208  -63.0571328  18.0508181                             Sint Maarten
# 209   55.4760200  -4.6609782                               Seychelles
# 210   38.5078794  35.0254725                                    Syria
# 211  -71.9738572  21.8304742                 Turks and Caicos Islands
# 212   18.6449321  15.3333334                                     Chad
# 213    0.9623212   8.5252980                                     Togo
# 214  101.0028835  15.1181635                                 Thailand
# 215   71.0136331  38.5304512                               Tajikistan
# 216   59.3709993  39.1155402                             Turkmenistan
# 217  125.8443864  -8.8288854                               East Timor
# 218 -174.8097544 -20.4281414                                    Tonga
# 219  -61.2656580  10.4573576                      Trinidad and Tobago
# 220    9.5528829  34.1195666                                  Tunisia
# 221   35.1689553  39.0616034                                   Turkey
# 222  120.9542726  23.7540034                                   Taiwan
# 223   34.8131021  -6.2756581              United Republic of Tanzania
# 224   32.3690754   1.2746902                                   Uganda
# 225   31.3832579  48.9965678                                  Ukraine
# 226  -56.0180680 -32.7995198                                  Uruguay
# 227 -112.4616706  45.6795520                 United States of America
# 228   63.1400180  41.7555408                               Uzbekistan
# 229   12.4338626  41.9017410                                  Vatican
# 230  -61.2012944  13.2247184         Saint Vincent and the Grenadines
# 231  -66.1818409   7.1242240                                Venezuela
# 232  -64.4714024  18.5259338                   British Virgin Islands
# 233  -64.8030050  17.9550140             United States Virgin Islands
# 234  106.2991386  16.6460145                                  Vietnam
# 235  167.6863674 -16.2263059                                  Vanuatu
# 236   35.2478213  31.9479803                                West Bank
# 237 -177.3481331 -13.8872499                        Wallis and Futuna
# 238 -172.1648148 -13.7532734                                    Samoa
# 239   47.5867601  15.9092765                                    Yemen
# 240   25.0838838 -29.0003377                             South Africa
# 241   27.7747508 -13.4582454                                   Zambia
# 242   29.8514503 -19.0042059                                 Zimbabwe
# 243  178.5198949  -7.7600583                                   Tuvalu
# 244  -53.2484462   3.9259572                            French Guiana