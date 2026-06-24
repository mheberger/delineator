"""
Constants for the delineator PYthon package.
"""
from importlib.resources import files

# Change this version number when the data files are updated. 
DATA_VERSION = "2026-06-22"   # bump ONLY when the data files change, independent of package version

# The width of half a pixel on the MERIT-Hydro grid, in decimal degrees
HALF_PIXEL = 0.000416667

#
PIXEL_AREA = 0.00083333 * 0.00083333

# MERIT-Hydro flow direction uses the ESRI standard for flow direction
DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)

# The scripts will try to use the spatialite extension if it is installed and fallback to geopandas otherwise.
# (This was mostly for testing during development.)
USE_SPATIALITE = False

VALID_MEGABASINS = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42,
                    43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74,
                    75, 76, 77, 78, 81, 82, 83, 84, 85, 86]

MEGABASINS_DB_FILE = files('delineator').joinpath('data', 'megabasins.db')

# List of *possible* output formats. Not all formats may be supported by the user's system
# and installed libraries.'
SUPPORTED_OUTPUT_FORMATS = ["shp", "gpkg", "geojson", "json", "kml", "parquet", "geoparquet"]

HASHES = {
    "vector/basins11.db": "81fdc13a635ab33faafe4b64e4296bb91dc3be84b88d50251010fd0ecfea3b35",
    "vector/basins12.db": "406454d4c28a487c21f47c44600c93be1458ca0b464f9b035c87490f7295b06e",
    "vector/basins13.db": "d27aa16b7bc8bbe968ed7d63d3e4525efdc0b24cf6cdd8a5589ebb8c35c1ecea",
    "vector/basins14.db": "19ac6ebcb1242fef40a7f29cbafcfb7fa7c785c29fc786e82fdee6d7d94356a3",
    "vector/basins15.db": "e031beb23f084ddb7c8feb3498e9f25b491f06e6bafffee3d5d4b49bc0968329",
    "vector/basins16.db": "7520e333ef11136f4c3e7c8101fe421ea308abc3ba7bc4ab2cc01086180aa34b",
    "vector/basins17.db": "22969f74031c32efb9873234ea74e0a678720cbf8b416c5d6db82d5bcef0fc86",
    "vector/basins18.db": "a0a9bdbbb2c010f4e6e6ae92d86c43a45b5cd2fe777329761693d6d8ff02a366",
    "vector/basins21.db": "7a8f7a11ecdc2610380c30994fed58907a1b2cf5088f769107d14daafe9781da",
    "vector/basins22.db": "f521a70d395b5f17b237151df9039f5748709cdf948cdc7d999d331c80c6bde3",
    "vector/basins23.db": "8c8737502e225939af697642274f7548c3e52869793289d60be686407f32f432",
    "vector/basins24.db": "da139e8552945e81e1d25dea0f9c5820d060e2ecf7c87ac28d0820d8ab4fd434",
    "vector/basins25.db": "29d86c1bf0fd5f263a5f948c7313b94b322f95b6c57a293294bf1b6072cf7fb3",
    "vector/basins26.db": "1974901b10c3ca62c637b91f24213634a42595387a92437e582d24fad873dd50",
    "vector/basins27.db": "b31666a1e1f7f68a361b9304d4ddc765eb20c8622f716233def7616dd9546838",
    "vector/basins28.db": "335eb6fee045eb7c735ce65459cec6013b7363136889adae96a3a5b31a9f5e50",
    "vector/basins29.db": "009b2f5de180d4c50a9f8b93dee6771048b9f4ce31b78f794f7976a2e1f36911",
    "vector/basins31.db": "95928c7663796a78c1745f6f35f36e8b8b91303c540e38749f276e42b2e91c69",
    "vector/basins32.db": "280e5fe4cfd67eae0628c86887877bcfaadbb75c0bd6ac960c915f7203d9764a",
    "vector/basins33.db": "3d23de4d013a3bc9eb2413c5389382878a76a1fe6f3f543a77440c624a40faea",
    "vector/basins34.db": "3cd01b60b45d6ea9faab2edf0c6eddfca2c3c418d4792b09da2ceecd366a8b75",
    "vector/basins35.db": "e171429b344b7e4fe928b9bc50281761849758b9d2c000d35b3cc8124a27ad41",
    "vector/basins36.db": "0756b196dfacfa4ecdb768f9362b183114bb6eb399e77fb16b19972ddd115d88",
    "vector/basins41.db": "eb8b31eb009f237f7237cf9c54bfc844f53953261c86c7a9cfb2fe8fb811420d",
    "vector/basins42.db": "1aceef20f908b7aafe26262252a24a5d3598d4f34727f3297216507a46354acb",
    "vector/basins43.db": "7d8ed6351dd0074024da20380f5d1f05f9d4bd7187017b5f30eb9e62704882ef",
    "vector/basins44.db": "edb900eebcf888f4b55bfde0d852da621c51d625bb4ca385ddbfe9ce725d02d0",
    "vector/basins45.db": "13dfd4fa896e26f46f405d448f5143289ecf9460ca7430463d6e65abb535de9a",
    "vector/basins46.db": "07796559518065f9331f3d2e61c079583673042862c93b883307c66fca67d69a",
    "vector/basins47.db": "f02cf47e8c55c1a0d1119466d4215bc4ce7d765f7fc8c533b49a5c3c16475362",
    "vector/basins48.db": "d57403a215dd50517e5237a092b711d1b1207416fce7b95b747386c8487f934d",
    "vector/basins49.db": "5d4315a0376a628ef92b74a4e199d77f5ab6b22626da3bba9f5f5eeba3ab5674",
    "vector/basins51.db": "db227655e91e2ef08ddcfa18ead5e7c4dda05c4050f336410d68350916bc65bc",
    "vector/basins52.db": "717af0a6fa88ca9705f133ae4e46e532dbd1121f56c8fdff7612741d7bdbfb80",
    "vector/basins53.db": "fc03933df4e4fa742b9fd68927834863d84481137ab64f46099adf04d258cf39",
    "vector/basins54.db": "f97eadddf095e2b8604f3fcc5fda4ca1646761fa88a6596fedaf5f4d54346c2a",
    "vector/basins55.db": "95580deae36af956269e3ee80874403e0feb9852f7f778d16d30cc298bae8e05",
    "vector/basins56.db": "c0330894ce917b63965310be5d2f2b412ee51a583febfff1030e165a81938dac",
    "vector/basins57.db": "8855d7021174c7a6925b64b25bae2c348634d2ea11508766069373d32b2a59b3",
    "vector/basins61.db": "0a5889d4655afac2a5fd9b4c6abdd992571b91cb996f2011f34074c41ea9e679",
    "vector/basins62.db": "8d58cb30d51d459a19470c8f92df9455a5ad7101273fc4f7cdafc71e6c956a73",
    "vector/basins63.db": "a2e91217cbc423cccf94c7a0d8f2b8d0397c3259d156b97eeabd532eeaa50497",
    "vector/basins64.db": "a623ed856a0f27a589a2aac33f14f2b44f4945842cbe5ca767a7c90a95b39800",
    "vector/basins65.db": "5ba6898363c1237ab44d7cd133139b02282f05307954e3f1394a2888cfcbeb51",
    "vector/basins66.db": "01f70c028d953e49d807fee1e5a22bcb70df29c4f308f2717120bdf954fcb611",
    "vector/basins67.db": "04cdca07c37997c809f2b742554ec3db285109ac07ea282182db2a395e4d156c",
    "vector/basins71.db": "6a3f9847f6aa41c08b364608a6e14a7914f1940ca7146a26c6e49965a46c78a4",
    "vector/basins72.db": "5c97fee6e979804f17d8f2f51f7ee4493bf522bfd0b24eb1e8c8960cc501a066",
    "vector/basins73.db": "7958655d8783d519aca3c0e6cabf0ecf4c078e65183c77ca83295a9e3e2dfc2e",
    "vector/basins74.db": "0969d1203dffc212d33a9ab386d6fe3d1b4274f2d1aeaddc47e8c9e6a3c9356b",
    "vector/basins75.db": "f1c4a9b5f2f8ddb82c0d7f6292a500e24839c0c0885ef2bcf159fb578aca8b5c",
    "vector/basins76.db": "60e107b73eead5fb84abfcf83c2b17ceee11fda754178218b05d84e2217d66f9",
    "vector/basins77.db": "903913221f2023c0affd2f13cf595f3c4fb8f1ddda77ed2361b4c3ea8607a77b",
    "vector/basins78.db": "c64cc3a35f2dd4b9e44f71275e3c9f5b0e2a5d5be17bf0bd5b672c54a22245b9",
    "vector/basins81.db": "11ac65a528beeb3b34549239c11e7c9fd6e55d80fa51ef90e2f6733c77fd5b88",
    "vector/basins82.db": "ef8e5290fa18b0097c6301ed812f905db238621308ed0eadbc0b69e93978dbab",
    "vector/basins83.db": "d01f326bc788df6c00ed8fb412720051c1196024d5b9fe470cef08e673c6d179",
    "vector/basins84.db": "e9dc564f6f2c1f88732cd0f721c3fd867bdd6e671044cee228460676d5759283",
    "vector/basins85.db": "859a8dabe6a467aa7a31f0ae6b333137d87a89f7e6309f99bc3ce605b6bfc4a7",
    "vector/basins86.db": "4b424de3c0fe40d4a455e766e9e9c9653a540b94ad6ba1a429e2198b204d7dc9",
    "vector/rivers11.db": "e395babeea01ca246f9aa98f5184b10cdb90a75889338ce89f1e87285b8780ed",
    "vector/rivers12.db": "288981d65a3c880071e751144ec1bfcf2cc0528cc65ba52f50e1f416f27c4e21",
    "vector/rivers13.db": "7f3dc00978a27096420d4cf1714cc940b941ec75994ab1dbf78e8e2b015d5f3a",
    "vector/rivers14.db": "5ae6be06798aa7f770ad91e3ca19e119f4f70e6a32f36e12b7463a4a70882c84",
    "vector/rivers15.db": "1f284ccb39e2c62134b3f33fef9abe9c0676fe67f78a767faab6f698dab1fced",
    "vector/rivers16.db": "5211be04c3ecfdb888d9794221a1075b464c3bfc3b3cdc3d9fc74e1a5ab35594",
    "vector/rivers17.db": "d3aa6fe8a074e39803ba1718453005cde3575a8231ed862f9c22e24b49a35eaf",
    "vector/rivers18.db": "36b51cd99bbec1f3bfc3abff5dedb92478dd41ef283d5c1ad4f63acb866ecaeb",
    "vector/rivers21.db": "b4d02dd21a47e72e26216cf2f4a58c3f72c4394d984d9ffc2d225d3a7f2e583c",
    "vector/rivers22.db": "5df3eb604fed2ec4e39415461703942ebe48f8adf6f1a817c05bd688e71d70ee",
    "vector/rivers23.db": "bd34bc3db35c274c9e8fb6e1f6fd1e6c9edaa5ea2db5dbd45d49de207e82b078",
    "vector/rivers24.db": "a40540d14b60507b3e9a49743346fab8dfe76d3d5be6e0bffa8e089e03eb1712",
    "vector/rivers25.db": "666f559fb4d78b43491313f9cc7a2cb1301f4ad446ba20989d05c006d8474521",
    "vector/rivers26.db": "e9310a6052f7e3a125d7acac974678d5027942fc637ed57fea5cd46fb20a52a0",
    "vector/rivers27.db": "2eb34cbfea5e7d3ce95566a5c107544f9a7342c6acd01271acb0bdbe69b38bfb",
    "vector/rivers28.db": "76315d7986698b483b816b6c56d26d44c51a34a21c1c73ac93fbcd5b1f9efa1f",
    "vector/rivers29.db": "8271023bc5f3e1bd7d9be5833bb6b7f795c3d68088f1ff3ec15e652e9032c3e9",
    "vector/rivers31.db": "eb1d0b82b8ffdf62970812840421386ec1775b3effdcbb86a7d45ff840c9ffb4",
    "vector/rivers32.db": "3f850d502a8a3a3e88287fc2fa2960bc15332caa038538475511287ca124fb0f",
    "vector/rivers33.db": "38f11ca2175600973d1ba9ffd01c5c14b30620651642b10d10bf5340f099a7bb",
    "vector/rivers34.db": "e3128b85d47d5e1c980ffe7377449a2a6b3caff0eec78e03b26cd0bf83b58589",
    "vector/rivers35.db": "9c83d711651431e8395402f27fc9c5f75a2395ff86a6a69e5e207629fb6fb5cb",
    "vector/rivers36.db": "2cc6c4a0bb97c90a60e7c44301b6f6839dfdf7c7971e3b23007d659117fa5109",
    "vector/rivers41.db": "910aa50867390ebb594f7f63141db847f4a47ca955fe4f6f2fe81d59748e0399",
    "vector/rivers42.db": "3ad08a5b007bd44a1ff064f50c5a3e502d1b1da6eafc0affd7af1951d35b49a2",
    "vector/rivers43.db": "3ea1570199e9974868613f48996e83848c140c087501c2d7d451ab31b3fde66c",
    "vector/rivers44.db": "62c639979501534d736e9b1e30f3de0ace43c0670d680d9201d1668b7bda4d1f",
    "vector/rivers45.db": "671ea0803b9180a2a7fbf7415f2aea3dfd451bc9d55b3188370b9b7c8156ee80",
    "vector/rivers46.db": "e91ef28c580d85b24cdbe68cf085d29560f5831752a4f7cc0817483381d62371",
    "vector/rivers47.db": "53431525adbe0e206571ec058bbc829750a6356086645ac80cd9542893ffbf58",
    "vector/rivers48.db": "959f60b31ca25f6e79b9c39cf188c982168150484a9e145e74cc4b508630538b",
    "vector/rivers49.db": "30471fed88912c8cdc0721d965ced432c6fec9e6d3776c05be5ec21b4254ada3",
    "vector/rivers51.db": "33b925cfb4c59859eede36d310b2c43612d436de12cbf8cd74ad0e2301790790",
    "vector/rivers52.db": "56188db23ca4044969e21030ee0602f5b27827bb5fa1b98f9ca607b92f4b1c43",
    "vector/rivers53.db": "fd226d5cbf117878517c4522e745d81bf7743319ea7ade01e6537ec4032bbd0d",
    "vector/rivers54.db": "c10a4e78bcc2247c11fbee67c22a5fc23d8ccd26047a7741030f23d131ac16ce",
    "vector/rivers55.db": "f8c700c777c9faef4b18aa9003b5f68fad0e445ca58985f058ac5e091404af73",
    "vector/rivers56.db": "0caed9b8c7bb877acda3daa9261b0ab269de7b7ce703052b6924e83f4857a151",
    "vector/rivers57.db": "fdff7e5e82357f64028d9b6b950d7e71c5f6b7b8bc98a02bb64f0741c20bfc83",
    "vector/rivers61.db": "30bccb97e9b08bae79d580b57eb78198861464e2662b0bee95bdae3b8b14ced7",
    "vector/rivers62.db": "863f1e13bea98f343765414f7c0d7442d8595480cd90c9f8ddfde67df56292c4",
    "vector/rivers63.db": "f41c8241a09c379d3d91cd044b0942ebf4079c153383221640d0d72c4191dd38",
    "vector/rivers64.db": "56f11316e505dbacbb998c990276b807a7ad2693c42ad0820243c28c62f881ae",
    "vector/rivers65.db": "701cb64ad7e3be7321bf34d9944b69f8bb4d9a527c23229b2b06dcd957d43745",
    "vector/rivers66.db": "99dc6ee86c419cfed55d13e58ecbf46d4bb4dfe763ee1bd387314569ddda4fd0",
    "vector/rivers67.db": "2d0a4d3411b6ef5c5eb0d1dbb1d9bb1d87fa8d6cd83d3f928708964127b5f962",
    "vector/rivers71.db": "d143d4923a2b6d014ee364498e6576bdf7d1d579bb6d35a2581033b729c976cb",
    "vector/rivers72.db": "2aeee034f3351b5e00d4982ea44d987bbb790679d48332039a58f6166a22da25",
    "vector/rivers73.db": "ca35a3599a5310e857ef57a28fbe2c6b143d09e75382f132d8bf158e5b0d51e8",
    "vector/rivers74.db": "cce66d3eccdd30bd16aa199b8439f08f2916eca80b88b09c5e7db8bdd48afee3",
    "vector/rivers75.db": "1aa0a9934ee81dcdfb6b4fc9d87c661af676e184962357dc5535408823823dcb",
    "vector/rivers76.db": "d8b465385359daae644a6bf935a2b77fe008dd48278dec50ddeea7f55ef95133",
    "vector/rivers77.db": "153848e68c5bbd6974061393539c28a8e152daaae5d00202d02ff3a1386a8f0e",
    "vector/rivers78.db": "effc6852743473273011cc584cc2e2de234fddff0ed023a75362e32414779ada",
    "vector/rivers81.db": "9c84e0a98e9ec059d40e5b9b9a630f1aabd8dd3f48d8419a406042cc28ac2320",
    "vector/rivers82.db": "1f3ec2645f30190b6c111570652ca5121d5c0028f123d6567b5dc689f9547239",
    "vector/rivers83.db": "c88c91be511e39655468dc83c66879edbab2e681753528311ae54194e4915f33",
    "vector/rivers84.db": "ae9ab6582e344bbabc446c52b0aa7c3e786a7249f54227536070bb7fc6ebacb1",
    "vector/rivers85.db": "704951c2abc6965ad88445b9576bfbfdc9d3574056a8cab5e3e836b28f5777d7",
    "vector/rivers86.db": "df14d0d979fa3ddcf25f61c75a742b6681f6dcd0560bc3400fe2a5a304141e02",
}
