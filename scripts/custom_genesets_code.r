################################################################################
### Custom Gene Sets from your curated file
### Add this code to your Shiny app
################################################################################

# In the UI section, replace the gene sets buttons with:

h5("Predefined Gene Sets:"),
fluidRow(
  column(12,
    h6("ER Stress & UPR:"),
    actionButton("load_upr_er", "UPR-ER", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_upr_mito", "UPR-Mito", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_erad", "ERAD", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_ire1_xbp1", "IRE1/XBP1s", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_perk", "PERK", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_atf6", "ATF6", class = "btn-info btn-sm", style = "margin: 2px;")
  )
),
br(),
fluidRow(
  column(12,
    h6("Cell Death & Survival:"),
    actionButton("load_cell_death", "Cell Death", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_autophagy", "Autophagy", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_prosurvival", "Pro-survival", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_apoptosis_classic", "Apoptosis", class = "btn-info btn-sm", style = "margin: 2px;")
  )
),
br(),
fluidRow(
  column(12,
    h6("Beta Cell Function:"),
    actionButton("load_beta_tfs", "Beta Cell TFs", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_insulin_secretion", "Insulin Secretion", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_ca_handling", "Ca2+/K+ Handling", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_disallowed", "Disallowed Genes", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_exocytosis", "Exocytosis", class = "btn-info btn-sm", style = "margin: 2px;")
  )
),
br(),
fluidRow(
  column(12,
    h6("Metabolism:"),
    actionButton("load_glycolysis", "Glycolysis/TCA", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_serine_1c", "Serine/1C", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_purine", "Purine Biosyn", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_redox", "Redox Stress", class = "btn-info btn-sm", style = "margin: 2px;")
  )
),
br(),
fluidRow(
  column(12,
    h6("Other Pathways:"),
    actionButton("load_ier", "Immediate Early", class = "btn-info btn-sm", style = "margin: 2px;"),
    actionButton("load_hyperins_mody", "Hyperins/MODY", class = "btn-info btn-sm", style = "margin: 2px;")
  )
)

# In the SERVER section, add these gene set handlers:

# UPR-ER genes (comprehensive)
observeEvent(input$load_upr_er, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Hspa5", "Atf6", "Atf6b", "Creb3l1", "Creb3l2", "Creb3l3", "Dnajc3",
      "Add1", "Apaf1", "Bak1", "Bcap31", "Dnajb2", "Dnajc10", "Amfr",
      "Ankzf1", "Aup1", "Bag6", "Brsk2", "Calr", "Calr3", "Calr4", "Canx",
      "Ccdc47", "Clgn", "Derl1", "Derl2", "Derl3", "Dnajb12", "Ecpas",
      "Edem1", "Edem2", "Edem3", "Erlec1", "Erlin1", "Erlin2", "Faf1",
      "Faf2", "Fbxo17", "Fbxo2", "Fbxo27", "Fbxo44", "Fbxo6", "Foxred2",
      "Get4", "Gmppb", "Herpud1", "Herpud2", "Hm13", "Jkamp", "Jkampl",
      "Man1a1", "Man1a2", "Man1b1", "Man1c1", "Marchf6", "Nccrp1", "Nploc4",
      "Nrros", "Os9", "Psmc6", "Rcn3", "Rhbdd1", "Rnf103", "Rnf121",
      "Rnf139", "Rnf185", "Rnf5", "Sdf2l1", "Sec61b", "Sec61bl", "Sel1l",
      "Sel1l2", "Selenos", "Sgta", "Stt3b", "Stub1", "Syvn1", "Tmem129",
      "Tmem67", "Tmub1", "Tmub2", "Tor1a", "Trim13", "Trim25", "Ube2g2",
      "Ube2j1", "Ube2j2", "Ube4a", "Ube4b", "Ubqln1", "Ubqln2", "Ubxn10",
      "Ubxn4", "Ubxn6", "Ubxn8", "Ufd1", "Uggt1", "Uggt2", "Umod", "Usp19",
      "Vcp", "Wfs1", "Yod1", "Dnajb9", "Ern1", "Ern2", "Map3k5", "Xbp1",
      "Eif2ak3", "Ddit3", "Atf4", "Eif2a", "Trib3", "Txnip", "Ppp1cc",
      "Ppp1r15a", "Ppp1r15b", "Manf", "Creld2"
    ), collapse = "\n"))
})

# UPR-Mito genes
observeEvent(input$load_upr_mito, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c("Atf5", "Lonp1", "Ubl5"), collapse = "\n"))
})

# ERAD genes
observeEvent(input$load_erad, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Hspa5", "Dnajb2", "Dnajc10", "Amfr", "Ankzf1", "Aup1", "Bag6",
      "Brsk2", "Calr", "Calr3", "Calr4", "Canx", "Ccdc47", "Clgn",
      "Derl1", "Derl2", "Derl3", "Dnajb12", "Ecpas", "Edem1", "Edem2",
      "Edem3", "Erlec1", "Erlin1", "Erlin2", "Faf1", "Faf2", "Fbxo17",
      "Fbxo2", "Fbxo27", "Fbxo44", "Fbxo6", "Foxred2", "Get4", "Gmppb",
      "Herpud1", "Herpud2", "Hm13", "Jkamp", "Jkampl", "Man1a1", "Man1a2",
      "Man1b1", "Man1c1", "Marchf6", "Nccrp1", "Nploc4", "Nrros", "Os9",
      "Psmc6", "Rcn3", "Rhbdd1", "Rnf103", "Rnf121", "Rnf139", "Rnf185",
      "Rnf5", "Sdf2l1", "Sec61b", "Sec61bl", "Sel1l", "Sel1l2", "Selenos",
      "Sgta", "Stt3b", "Stub1", "Syvn1", "Tmem129", "Tmem67", "Tmub1",
      "Tmub2", "Tor1a", "Trim13", "Trim25", "Ube2g2", "Ube2j1", "Ube2j2",
      "Ube4a", "Ube4b", "Ubqln1", "Ubqln2", "Ubxn10", "Ubxn4", "Ubxn6",
      "Ubxn8", "Ufd1", "Uggt1", "Uggt2", "Umod", "Usp19", "Vcp", "Wfs1",
      "Yod1", "Dnajb9", "Eif2ak3"
    ), collapse = "\n"))
})

# IRE1/XBP1s pathway
observeEvent(input$load_ire1_xbp1, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Jun", "Traf2", "Dnajb9", "Ern1", "Ern2", "Map3k5", "Xbp1"
    ), collapse = "\n"))
})

# PERK pathway
observeEvent(input$load_perk, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Eif2ak3", "Ddit3", "Atf4", "Eif2a", "Trib3", "Txnip",
      "Ppp1cc", "Ppp1r15a", "Ppp1r15b"
    ), collapse = "\n"))
})

# ATF6 pathway
observeEvent(input$load_atf6, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Hspa5", "Atf6", "Atf6b", "Creb3l1", "Creb3l2", "Creb3l3", "Dnajc3"
    ), collapse = "\n"))
})

# Cell Death genes (comprehensive)
observeEvent(input$load_cell_death, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Mapk8", "Pmaip1", "Traf2", "Bcl2l11", "Bid", "Casp3", "Casp9",
      "Dnm1l", "Lmna", "Tnfrsf10b", "Acin1", "Apc", "Apip", "Aven",
      "Bad", "Bax", "Bcl2l1", "Birc2", "Bmf", "Bmx", "Casp6", "Casp7",
      "Casp8", "Cd14", "Cflar", "Clspn", "Cycs", "Dcc", "Dffa", "Dffb",
      "Diablo", "Dsg1a", "Dsg2", "Dsg3", "Dsp", "Dynll1", "Dynll2",
      "Fadd", "Fas", "Fasl", "Fnta", "Gas2", "Gm10053", "Gsdmd", "Gsdme",
      "Gsn", "Gzmb", "H1f0", "H1f1", "H1f2", "H1f4", "H1f5", "Hmgb1",
      "Hmgb2", "Kpna1", "Kpnb1", "Lmnb1", "Ly96", "Mapk1", "Mapk3",
      "Mapt", "Nmt1", "Ocln", "Oma1", "Opa1", "Pkp1", "Plec", "Ppp3cc",
      "Ppp3r1", "Prkcd", "Prkcq", "Ptk2", "Ripk1", "Rock1", "Satb1",
      "Septin4", "Sfn", "Sorbs2", "Sptan1", "Stk24", "Stk26", "Ticam1",
      "Ticam2", "Tjp1", "Tjp2", "Tlr4", "Tnfsf10", "Tradd", "Vim",
      "Xiap", "Ywhab", "Ywhae", "Ywhag", "Ywhah", "Ywhaq", "Ywhaz",
      "Ctnnb1", "Add1", "Apaf1", "Bak1", "Bcap31"
    ), collapse = "\n"))
})

# Autophagy genes
observeEvent(input$load_autophagy, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Psen1", "Atg12", "Atg13", "Atg14", "Atg16l1", "Atg16l2", "Atg2a",
      "Atg3", "Atg5", "Atg7", "Atp6ap2", "Becn1", "Calcoco2", "Epg5",
      "Lamp1", "Pik3c3", "Rb1cc1", "Sqstm1", "Ulk1", "Ulk2", "Uvrag"
    ), collapse = "\n"))
})

# Pro-survival genes
observeEvent(input$load_prosurvival, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Psen1", "Dad1", "H13", "Plk", "Rfc4", "Ssr4", "Manf"
    ), collapse = "\n"))
})

# Classic apoptosis
observeEvent(input$load_apoptosis_classic, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Casp3", "Casp7", "Casp8", "Casp9", "Bax", "Bcl2", "Bcl2l1",
      "Mcl1", "Parp1", "Cycs", "Bid", "Bad", "Bak1", "Apaf1"
    ), collapse = "\n"))
})

# Beta cell transcription factors
observeEvent(input$load_beta_tfs, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Ctnnb1", "Foxa2", "Neurod1", "Nsd1", "Pdx1", "Bmp5", "Crtc2",
      "Ep300", "Foxo1", "Ins1m", "Isl1", "Mafa", "Mafb", "Mnx1",
      "Nkx2-2", "Nkx6-1", "Pax6", "Rfx2", "Rfx3", "Rfx6", "Rreb1",
      "Scrt1", "Tcf7l2", "Tshz1", "Ucn3", "Wnt4", "Zfp148"
    ), collapse = "\n"))
})

# Insulin secretion genes
observeEvent(input$load_insulin_secretion, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Scgn", "Vgf", "Doc2b", "Napa", "Nsf", "Stx1a", "Stx4", "Stxbp1",
      "Stxbp3", "Stxbp5", "Stxbp5l", "Syt4", "Syt7", "Syt9", "Unc13a",
      "Unc13b", "Vamp2", "Ap3b1", "Cdc42", "Chga", "Chgb", "Sec23a",
      "Sec24a", "Sec24d", "Hsp90b1", "Ptbp1", "Pam", "Pcsk1", "Pcsk2",
      "Scg5", "Rab27a", "Rab3a", "P4hb", "Ssr1", "Cdkal1", "Cpe",
      "Ddx1", "Eif3a", "Eif4b", "Erc1", "Ero1b", "Ero1lb", "Ins2",
      "Slc30a5", "Slc30a8", "Yipf5", "Abcc8", "Ano1", "Cacna1d",
      "Kcnh2", "Kcnq1", "Appl1", "Pgm1", "Gck", "Glud1", "Hadh",
      "Adcy7", "Trmt10a", "Pmm2", "Hnf1a", "Hnf1b", "Hnf4a", "Klf11",
      "Pax4", "Yars", "Blk", "Cel", "Eif2s3", "Glis3", "Ins", "Insr",
      "Ucp2", "Kcnj11", "Foxa2", "Neurod1", "Nsd1", "Pdx1"
    ), collapse = "\n"))
})

# Calcium and potassium handling
observeEvent(input$load_ca_handling, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Kcnj11", "Atp2a1", "Atp2a3", "Itpr1", "Atp2a2", "Atp2b1",
      "Atp2b2", "Atp2b3", "Atp2b4", "Atp2c1", "Cacna1c", "Itpr2",
      "Itpr3", "Letm1", "Mcoln1", "Mcoln2", "Mcoln3", "Mcu", "Mcub",
      "Micu1", "Micu2", "Micu3", "Orai1", "Orai2", "Orai3", "Pkd2",
      "Pkd2l1", "Ryr1", "Ryr2", "Ryr3", "Slc24a1", "Slc24a2", "Slc24a3",
      "Slc24a4", "Slc24a5", "Slc8a1", "Slc8a2", "Slc8a3", "Slc8b1",
      "Smdt1", "Trpc1", "Trpc2", "Trpc3", "Trpc4", "Trpc5", "Trpc6",
      "Trpc7", "Trpv1", "Trpv2", "Trpv3", "Trpv4", "Trpv5", "Trpv6",
      "Abcc8", "Ano1", "Cacna1d", "Kcnh2", "Kcnq1"
    ), collapse = "\n"))
})

# Disallowed beta cell genes
observeEvent(input$load_disallowed, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Hk1", "Slc16a1", "Acot7", "Arhgdib", "Cat", "Cxcl12", "Hsd11b1",
      "Igf1", "Igfbp4", "Itih5", "Ldha", "Maf", "Mgll", "Oat", "Pdgfra",
      "Smad3", "Smoc2", "Yap1", "Zfp36l1", "Zyx"
    ), collapse = "\n"))
})

# Exocytosis genes
observeEvent(input$load_exocytosis, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Scgn", "Vgf", "Doc2b", "Napa", "Nsf", "Stx1a", "Stx4", "Stxbp1",
      "Stxbp3", "Stxbp5", "Stxbp5l", "Syt4", "Syt7", "Syt9", "Unc13a",
      "Unc13b", "Vamp2", "Ap3b1", "Cdc42", "Chga", "Chgb"
    ), collapse = "\n"))
})

# Glycolysis and TCA cycle
observeEvent(input$load_glycolysis, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Aldoa", "Aldob", "Aldoc", "Bpgm", "Eno1", "Eno2", "Eno3",
      "Gapdh", "Gapdhs", "Gpi", "Hk2", "Hkdc1", "Pfkfb1", "Pfkfb2",
      "Pfkfb3", "Pfkfb4", "Pfkl", "Pfkm", "Pfkp", "Pgam1", "Pgam2",
      "Pgam4", "Pgk1", "Pklr", "Pkm", "Ppp2ca", "Ppp2r1a", "Ppp2r1b",
      "Ppp2r5d", "Tpi1", "Aco2", "Acsl3", "Acsl4", "Adsl1", "Cs",
      "Dld", "Dlst", "Fh", "Glud2", "Gucy2c", "Idh1", "Idh2", "Idh3a",
      "Idh3b", "Idh3g", "Mdh2", "Me3", "Nnt", "Ogdh", "Sdha", "Sdhb",
      "Sdhc", "Sdhd", "Sucla2", "Suclg1", "Suclg2", "Gck", "Glud1",
      "Hadh", "Ppp2cb"
    ), collapse = "\n"))
})

# Serine/1C metabolism
observeEvent(input$load_serine_1c, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Aldh1l2", "Mthfd1", "Mthfd2", "Phgdh", "Psat1", "Psph",
      "Shmt1", "Shmt2"
    ), collapse = "\n"))
})

# Purine biosynthesis
observeEvent(input$load_purine, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Adsl", "Adss", "Ak1", "Ampd2", "Ampd3", "Aprt", "Atic", "Gmps"
    ), collapse = "\n"))
})

# Redox stress
observeEvent(input$load_redox, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Bcl2", "Dusp10", "Elk1", "Map2k3", "Map2k4", "Map2k6", "Mapk11",
      "Mapk12", "Mapk13", "Mapk14", "Mapk9", "Max", "Mef2c", "Mknk1",
      "Mknk2", "Myc", "Pla2g4a", "Stat1", "Txn", "Mapk8", "Jun", "Ddit3"
    ), collapse = "\n"))
})

# Immediate early response
observeEvent(input$load_ier, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Jun", "Atf3", "Ccl2", "Bhlhe40", "Ccn1", "Ccn2", "Ccnl1",
      "Cebpd", "Csrnp1", "Cxcl3", "Dusp1", "Dusp5", "Dusp6", "Egr1",
      "Egr3", "F3", "Flg", "Fos", "Fosb", "Gadd45b", "Gbp1", "Gem",
      "Hbegf", "Ier3", "Il6", "Junb", "Klf10", "Klf6", "Klhl21",
      "Ldlr", "Mcl1", "Nfkbia", "Nfkbiz", "Nr2c2", "Nr4a1", "Nr4a2",
      "Plau", "Rcan1", "Sgk1", "Slc2a3", "Srf", "Tnfaip3", "Trib1",
      "Tsc22d1", "Zfp36", "Pmaip1"
    ), collapse = "\n"))
})

# Hyperinsulinism and MODY genes
observeEvent(input$load_hyperins_mody, {
  updateTextAreaInput(session, "gene_list",
    value = paste(c(
      "Kcnj11", "Foxa2", "Neurod1", "Nsd1", "Pdx1", "Hk1", "Slc16a1",
      "Abcc8", "Ano1", "Cacna1d", "Kcnh2", "Kcnq1", "Appl1", "Pgm1",
      "Gck", "Glud1", "Hadh", "Adcy7", "Trmt10a", "Pmm2", "Hnf1a",
      "Hnf1b", "Hnf4a", "Klf11", "Pax4", "Yars", "Blk", "Cel",
      "Eif2s3", "Glis3", "Ins", "Insr", "Ucp2"
    ), collapse = "\n"))
})