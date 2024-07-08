import matplotlib.pyplot as plt

def get_protein_length(idr_positions):
    return [(idr[4], idr) for idr in idr_positions]

def plot_phosphorylation_sites(idr_positions, title):
    # Sort the IDR positions by protein length in descending order
    idr_positions_sorted = sorted(get_protein_length(idr_positions), key=lambda x: x[0], reverse=False)
    idr_positions_sorted = [idr[1] for idr in idr_positions_sorted]

    fig, ax = plt.subplots(figsize=(10, 8))

    y_labels = []
    y_ticks = []

    for i, (start, end, protein, psite, protein_length) in enumerate(idr_positions_sorted):
        y = len(idr_positions_sorted) - i
        y_labels.append(protein)
        y_ticks.append(y)

        # Plot the IDR region
        centered_start = start - psite
        centered_end = end - psite
        ax.plot([centered_start, centered_end], [y, y], color='black', lw=1.5)
        centered_psite = 0

        ax.plot([centered_start, centered_end], [y, y], color='black', lw=1.5)
        ax.plot(centered_psite, y, 'ro', markersize=1.5)  # Red dot for the p-site

        # Plot start and end markers for the IDR
        ax.plot(centered_start, y, 'go', markersize=1.5)  # Green dot for start
        ax.plot(centered_end, y, 'bo', markersize=1.5)  # Blue dot for end

        # Plot protein length
        centered_length_start = -start + centered_start
        centered_length_end = centered_end + (protein_length - end)
        ax.plot([centered_length_start, centered_length_end], [y, y], color='gray', lw=1.5, linestyle='dotted')

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, rotation=0, ha='right', fontsize=5)
    ax.tick_params(axis='y', labelsize=3)
    ax.set_xlabel('Position')
    ax.set_title(title)
    plt.gca().invert_yaxis()  # To have the first protein on top

    # Create legends
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color='black', lw=3, label='IDR Region'),
        Line2D([0], [0], color='gray', lw=3, label='Total Protein Length', linestyle='dotted'),
        Line2D([0], [0], color='red', marker='o', linestyle='None', label='p-site'),
        Line2D([0], [0], color='green', marker='o', linestyle='None', label='IDR Start'),
        Line2D([0], [0], color='blue', marker='o', linestyle='None', label='IDR End')
    ]
    ax.legend(handles=custom_lines, loc='lower right')

    plt.show()


def main():

    idr_positions = [
        (1796, 2086, 'RHG32', 1796, 2087),
        # (1036, 1036, 'MYO1C', 736),
        # (1, 12, 'HSPB6', 16),
        (265, 330, 'CHK1', 296, 476),
        (284, 930, 'MYPT1', 472, 1030),
        # (121, 180, 'MDM4', 367),
        # (1, 87, 'PDPK1', 354),
        (83, 376, 'MEFV', 242, 781),
        (1, 117, 'B2L11', 87, 198),
        (224, 593, 'FOXO3', 253, 673),
        (787, 1103, 'KIF1C', 1092, 1103 ),
        (228, 354, 'TBCD4', 341, 1298),
        # (432, 505, 'MITF', 280),
        (1, 188, 'SNN', 44, 88),
        (1619, 2440, 'NCOR1', 2335, 2440),
        (274, 820, 'ULK1', 556, 1050 ),
        # (139, 259, 'GRAP2', 262),
        (353, 478, 'BEST1', 358, 585),
        (87, 106, 'SASH1', 90,1247 ),
        (67, 143, 'CDK14', 119, 469),
        (1, 35, 'SNAI1', 11, 264),
        (888, 1035, 'M3K6', 964, 1288),
        (1, 306, 'PAK4', 181, 591),
        (498, 1003, 'ABL1', 735, 1130 ),
        # (63, 177, 'ESR1'),
        (201, 341, 'RAF1', 259, 648),
        (87, 188, 'GCR', 134, 777),
        (1, 44, 'HG2A', 8, 296),
        (277, 393, 'P53', 366, 393 ),
        (728, 768, 'ITB2', 758, 769 ),
        (756, 798, 'ITB1', 788, 798),
        (1, 112, 'K1C18', 34, 430),
        (1, 96, 'TY3H', 71, 528 ),
        # (318, 484, 'GP1BA'),
        (1283, 1367, 'IGF1R', 1313, 1367),
        (323, 1106, 'GLI1', 640, 1106),
        (471, 565, 'WEE2', 557, 567),
        (535, 1538, 'GLI2', 640, 1586),
        (826, 1405, 'GLI3', 1006, 1580),
        (123, 297, 'ARAF', 214, 606),
        (1, 630, 'TAU', 531, 758),
        (771, 821, 'FGFR1', 779, 822),
        (1002, 1032, 'ITA4', 1011, 1032),
        (1, 50, 'GFAP', 8, 432),
        (295, 462, 'BRAF', 365, 766),
        # (92, 93, 'NCK1'),
        (318, 561, 'PRLR', 415, 622),
        # (105, 162, 'ZNF28'),
        (273, 382, 'CXA1', 373, 382),
        (1, 274, 'NELFE', 115, 380),
        # (1, 165, 'TFEB'),
        (702, 815, 'SL9A1', 703, 815),
        (238, 403, 'HNF1A', 247, 631),
        (1, 77, 'PHOS', 54, 246),
        (1, 427, 'CD11B', 113, 795),
        (773, 821, 'FGFR2', 782, 821),
        # (457, 475, 'TGM2'),
        (340, 360, 'NR4A1', 351, 598  ),
        # (1, 1, 'COF1'),
        (358, 490, 'PTN3', 359, 913),
        (1, 107, 'TTP', 60, 326),
        (48, 124, 'P85A', 83, 724),
        # (1, 284, 'WEE1'),
        (1, 330, 'MPIP1', 178, 524),
        (1, 387, 'MPIP2', 151, 580),
        (102, 228, 'MPIP3', 216, 473 ),
        (1100, 1229, 'L1CAM', 1181,1257 ),
        (411, 896, 'IL3RB', 601, 897),
        (551, 554, 'CTNB1', 552, 781),
        (241, 1242, 'IRS1', 374, 1242),
        (315, 765, 'UBP8', 718, 1118),
        (861, 1078, 'CASR', 895, 1078),
        (230, 270, 'AQP2', 256, 271 ),
        # (311, 361, 'CASP2'),
        (264, 359, 'KSYK', 297, 635),
        (352, 627, 'ACHA4', 467, 627),
        (1, 198, 'CDN1B', 157, 198),
        (1, 466, 'YAP1', 127, 504 ),
        # (921, 1157, 'NRIP1', ),
        # (272, 337, 'KC1A'),
        (1, 58, 'CENPA', 7, 140),
        (285, 363, 'PCY1A', 288, 367),
        # (1, 1, 'STAR'),
        (197, 646, 'NUMB', 295, 651),
        (654, 1067, 'RGS3', 943, 1198),
        (918, 977, 'TSC2', 939, 1807),
        (91, 271, 'ETV1', 216, 477),
        (178, 226, 'RND2', 223, 227),
        (327, 398, 'PLK1', 330, 603),
        (715, 812, 'ATX1', 775, 815),
        (1, 69, 'GEM', 23, 296 ),
        (87, 667, 'HDAC4', 350, 1084 ),
        # (65, 478, 'BCAR1'),
        # (428, 519, 'SIK1'),
        (2, 268, 'MLF1', 34, 268),
        # (1, 15, 'RND3'),
        (1, 47, 'H31', 29, 136),
        (133, 462, '3BP2', 278,561 ),
        (230, 508, 'SRPK2', 492, 688 ),
        # (653, 667, 'GPSM2'), no psite in IDR
        (173, 295, 'FOXO4', 197, 505),
        # (1, 116, 'CDK16'), no psite in IDR
        (99, 286, 'MDM2', 166, 491),
        (291, 402, 'KPCE', 346, 737),
        (765, 960, 'KIF23', 814, 960),
        (1359, 1400, 'RON', 1394, 1400),
        (48, 56, 'BTK', 51, 659),
        (44, 121, 'TISB', 92, 338 ),
        (402, 573, 'KLC1', 545, 573 ),
        (1015, 1333, 'SOS1', 1161, 1333),
        (1, 104, 'FOXO1', 24, 655),

        (1, 429, 'TIAM1', 172, 1591),
        (651, 653, 'TRI32', 652, 653),
        (74, 421, 'LCP2', 376, 533),
        (113, 327, 'ZBT17', 291, 803),
        (375, 828, 'MK07', 486, 816),
        (1, 94, 'SKP2', 72, 424),
        # (1, 169, 'GRB10'), no psite in IDR
        # (371, 439, 'CTBP1'), no psite in IDr
        (405, 601, 'TNK1', 502, 666),
        (211, 381, 'RIN1', 351, 783),
        (1141, 1323, 'SPTA2', 1302, 2472),
        (1, 481, 'SRC8', 298, 550),
        (400, 652, 'PDE3A', 428,1141 ),
        # (1, 123, 'BECN1'), no psite in IDR
        # (911, 1036, 'KANK1'), no psite in IDR
        (106, 459, 'LPIN1', 287, 890 ),
        (1, 378, 'NFAC4', 272, 902),
        # (911, 1036, 'NMDE3'), no IDR in psite
        (347, 494, 'RHG19', 422, 494),
        # (1, 44, 'KS6A1'), no psite in any IDR
        (314, 607, 'TESK1', 437, 626),
        (302, 433, 'STK11', 336, 433),
        (1, 48, 'SNAT', 31, 207),
        (1, 47, 'H31T', 29, 136),
        (124, 204, 'TPD53', 137, 204),
        # (346, 1133, 'RHG31'), no psite in IDR
        (668, 888, 'RGPA2', 715, 1873),
        (1, 626, 'CRTC2', 171, 693 ),
        (550, 693, 'SMAG2', 642, 694 ),
        # (384, 504, 'LRRK2'), no IDR

        (647, 1048, 'DAB2P', 728, 1189 ),
        (1, 195, 'TET2', 99, 2002),
        (371, 475, 'KLC3', 431, 504),
        (1, 219, 'RHDF2', 113, 856),
        (517, 625, 'LRFN4', 585, 635),

        (1093, 1286, 'RICTR', 1135, 1708),
        # (350, 642, 'PADI6'), no IDR
        (353, 711, 'MARK2', 596, 788),
        (625, 838, 'ATG9A', 761, 839),
        (237, 255, 'KCNKI', 252, 384),

        (1, 32, 'RHG22', 16, 698),
        (512, 724, 'TBCD1', 596, 1168),
        (144, 327, 'NOXA1', 172, 476 ),
        (173, 590, 'RIMS1', 357, 1692),
        (174, 469, 'PACS2', 437, 889),
        (280, 649, 'LSR', 493, 649),
        (172, 324, 'KSR1', 309, 923),
        # (1, 63, 'TPH2'), no psite in IDR
        (1, 122, 'REM2', 68, 340),
        # (838, 962, 'RPTOR'), no psite in IDr
        (1, 100, 'KKCC1', 74, 505),
        (440, 605, 'HJURP', 486, 748),
        (1, 115, 'ZNRF1', 50, 227),
        # (1, 60, 'CCNY'), no psite in IDR
        (1, 146, 'ZNRF2', 82, 242),
        (1, 535, 'HDAC7', 358, 952),

        (427, 932, 'SSH1', 834, 1049),
        (202, 447, 'RUBIC', 248, 972 ),
        # (150, 248, 'RND1'),no IDR
        (1, 168, 'BAD', 75, 168),
        (645, 986, 'ARHG2', 886, 986),
        (1, 256, 'AKTS1', 183, 256),
        # (62, 243, 'DDT4L'), no IDR
        (67, 230, 'EDC3', 161, 508),
        (1, 350, 'JUB', 230, 538 ),
        (124, 594, 'NED4L', 342, 975),
        (693, 696, 'SLIK1', 695, 696),
        (1, 214, 'CIC', 173, 2517),
        (27, 56, 'RMD3', 46, 470),
        # (1172, 1249, 'M3K5'),
        (128, 334, 'M3K3', 294, 626),
        (80, 346, 'PKP2', 82, 881),
        (60, 251, 'CDCA7', 163, 371),
        (287, 728, 'PP12C', 452, 782),
        (294, 335, 'NDEL1', 336, 345),
        (1, 359, 'WWTR1', 89, 400),
        (125, 240, 'MFF', 157, 342),
        (395, 619, 'KLC2', 582, 622),
        # (436, 602, 'SIK2'), no psite in IDr
        # (157, 167, 'ISCU'), no psite in IDR
        (1, 250, 'MYLK2', 151, 596),
        # (161, 251, 'REEP4'), no psite in IDR
        (317, 480, 'ZN395', 396, 513),
        (1065, 1162, 'CENPJ', 1109, 1338),
        # (262, 659, 'WBS14'), no psite in IDR
        # (1, 11, 'NGB'), no psite in IDr
        (87, 376, 'PAK6', 113, 681),
        # (182, 203, 'RASF1'), no psite in IDr
        # (46, 78, 'RGS18'), no psite in IDR
        (401, 619, 'KLC4', 554, 619),
        # (1, 72, 'DDIT4'), no psite in idr
        # (171, 172, 'DACT1'), no psite in idr
        (397, 550, 'DTL', 464, 730),
        (746, 900, 'ADA22', 857, 906),
        # (179, 654, 'TBCD7'), no idr
        (224, 320, 'PI4KB', 294, 816),
        (220, 381, 'ERRFI', 251, 462),
        (137, 373, 'DBNL', 291, 430),
        (253, 356, 'ING1', 342, 422),
        (263, 332, 'PITC1', 299, 332),
        (72, 650, 'HDAC9', 451, 1011),
        # (1, 1, 'PRP19'), no idr
        (456, 1091, 'SYNP2', 611, 1093),
        (2238, 3459, 'BSN', 2851, 3926),
        (165, 669, 'RIMS2', 366, 1411),
        (507, 846, 'EXO1', 746, 846),
        (255, 372, 'BAIP2', 340, 552),
        (239, 634, 'GAB2', 391, 676),
        (455, 683, 'HDAC5', 498, 1122),
        (79, 430, 'SH2B3', 150, 575),
        (142, 403, 'RP3A', 272, 694),
        (428, 1319, 'SIK3', 493, 1321),
        (134, 353, 'M3K2', 282, 525),
        # RNF 11 no psite in IDR
        (1, 72, 'CBY1', 20, 126),
        (297, 1338, 'IRS2', 1148, 1338),
        (1, 159, 'GPSM3', 35, 159),
        (80, 142, 'RN115', 133, 304),
        (1913, 2219, 'CABIN', 2126, 2219),
        (194, 609, 'NUMBL', 324, 609),
        (24, 435, 'LIMD1', 316, 676)







]

    plot_phosphorylation_sites(idr_positions, "Phosphorylation Sites Plotted on IDRs of 14-3-3 Client Proteins")



main()