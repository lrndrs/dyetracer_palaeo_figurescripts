import pandas as pd

DYE_TABLE = pd.DataFrame([
    # ------------------------------------------------------------
    # EIS (European Ice Sheet Complex)
    # ------------------------------------------------------------
    {
        "dye": "dye00",
        "region": "EIS_MedSea",
        "ice_sheet": "European Ice Sheet Complex",
        "volume_m3": 1.79e13,
        "group": "EIS",
    },
    {
        "dye": "dye01",
        "region": "EIS_BayOfBiscay",
        "ice_sheet": "European Ice Sheet Complex",
        "volume_m3": 4.92e12,
        "group": "EIS",
    },
    {
        "dye": "dye02",
        "region": "EIS_NorwegianSea",
        "ice_sheet": "European Ice Sheet Complex",
        "volume_m3": 5.85e12,
        "group": "EIS",
    },
    {
        "dye": "dye03",
        "region": "EIS_Arctic",
        "ice_sheet": "European Ice Sheet Complex",
        "volume_m3": 7.02e12,
        "group": "EIS",
    },

    # ------------------------------------------------------------
    # LAU (Laurentide Ice Sheet)
    # ------------------------------------------------------------
    {
        "dye": "dye04",
        "region": "LAU_Arctic",
        "ice_sheet": "Laurentide Ice Sheet",
        "volume_m3": 1.34e13,
        "group": "LAU",
    },
    
    
    # ------------------------------------------------------------
    # GIS (Greenland Ice Sheet)
    # ------------------------------------------------------------
    {
        "dye": "dye05",
        "region": "GIS_GreenlandSea",
        "ice_sheet": "Greenland Ice Sheet",
        "volume_m3": 4.49e12,
        "group": "GIS",
    },


    # ------------------------------------------------------------
    # LIS South 
    # ------------------------------------------------------------


    {
        "dye": "dye06",
        "region": "LAU_LabradorSea",
        "ice_sheet": "Laurentide Ice Sheet",
        "volume_m3": 1.33e13,
        "group": "LAU",
    },
    
    
    {
        "dye": "dye07",
        "region": "LAU_StLawrence",
        "ice_sheet": "Laurentide Ice Sheet",
        "volume_m3": 1.58e13,
        "group": "LAU",
    },
    {
        "dye": "dye08",
        "region": "LAU_GulfOfMexico",
        "ice_sheet": "Laurentide Ice Sheet",
        "volume_m3": 1.19e13,
        "group": "LAU",
    },

   
])