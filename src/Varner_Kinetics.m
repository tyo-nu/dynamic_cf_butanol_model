function [dx, rate_vector] = Kinetics(t,x, model)
    M_g6p_c = x(1);
	M_f6p_c = x(2);
	M_fdp_c = x(3);
	M_dhap_c = x(4);
	M_g3p_c = x(5);
	M_13dpg_c = x(6);
	M_3pg_c = x(7);
	M_2pg_c = x(8);
	M_oaa_c = x(9);
	M_coa_c = x(10);
	M_accoa_c = x(11);
	M_6pgl_c = x(12);
	M_6pgc_c = x(13);
	M_ru5p_D_c = x(14);
	M_xu5p_D_c = x(15);
	M_r5p_c = x(16);
	M_s7p_c = x(17);
	M_e4p_c = x(18);
	M_2ddg6p_c = x(19);
	M_cit_c = x(20);
	M_icit_c = x(21);
	M_akg_c = x(22);
	M_succoa_c = x(23);
	M_q8_c = x(24);
	M_fum_c = x(25);
	M_q8h2_c = x(26);
	M_mql8_c = x(27);
	M_mqn8_c = x(28);
	M_h_e = x(29);
	M_ppi_c = x(30);
	M_glx_c = x(31);
	M_actp_c = x(32);
	M_etoh_c = x(33);
	M_for_c = x(34);
	M_nh3_c = x(35);
	M_arg_L_c = x(36);
	M_h2s_c = x(37);
	M_thf_c = x(38);
	M_mlthf_c = x(39);
	M_aicar_c = x(40);
	M_5mthf_c = x(41);
	M_chor_c = x(42);
	M_h2o2_c = x(43);
	M_mglx_c = x(44);
	M_prop_c = x(45);
	M_indole_c = x(46);
	M_cadav_c = x(47);
	M_gaba_c = x(48);
	M_glycoA_c = x(49);
	M_78mdp_c = x(50);
	M_4adochor_c = x(51);
	M_4abz_c = x(52);
	M_78dhf_c = x(53);
	M_dhf_c = x(54);
	M_methf_c = x(55);
	M_10fthf_c = x(56);
	M_prpp_c = x(57);
	M_hco3_c = x(58);
	M_clasp_c = x(59);
	M_or_c = x(60);
	M_omp_c = x(61);
	M_5pbdra = x(62);
	M_gar_c = x(63);
	M_fgar_c = x(64);
	M_fgam_c = x(65);
	M_air_c = x(66);
	M_cair_c = x(67);
	M_saicar_c = x(68);
	M_faicar_c = x(69);
	M_imp_c = x(70);
	M_xmp_c = x(71);
	GENE_CAT = x(72);
	RNAP = x(73);
	OPEN_GENE_CAT = x(74);
	mRNA_CAT = x(75);
	RIBOSOME = x(76);
	RIBOSOME_START_CAT = x(77);
	M_ala_L_c_tRNA = x(78);
	M_arg_L_c_tRNA = x(79);
	M_asn_L_c_tRNA = x(80);
	M_asp_L_c_tRNA = x(81);
	M_cys_L_c_tRNA = x(82);
	M_glu_L_c_tRNA = x(83);
	M_gln_L_c_tRNA = x(84);
	M_gly_L_c_tRNA = x(85);
	M_his_L_c_tRNA = x(86);
	M_ile_L_c_tRNA = x(87);
	M_leu_L_c_tRNA = x(88);
	M_lys_L_c_tRNA = x(89);
	M_met_L_c_tRNA = x(90);
	M_phe_L_c_tRNA = x(91);
	M_pro_L_c_tRNA = x(92);
	M_ser_L_c_tRNA = x(93);
	M_thr_L_c_tRNA = x(94);
	M_trp_L_c_tRNA = x(95);
	M_tyr_L_c_tRNA = x(96);
	M_val_L_c_tRNA = x(97);
	PROTEIN_CAT = x(98);
	tRNA = x(99);
	M_glc_D_c = x(100);
	M_pep_c = x(101);
	M_pyr_c = x(102);
	M_ac_c = x(103);
	M_lac_D_c = x(104);
	M_mal_L_c = x(105);
	M_atp_c = x(106);
	M_adp_c = x(107);
	M_amp_c = x(108);
	M_gtp_c = x(109);
	M_gdp_c = x(110);
	M_gmp_c = x(111);
	M_utp_c = x(112);
	M_udp_c = x(113);
	M_ump_c = x(114);
	M_ctp_c = x(115);
	M_cdp_c = x(116);
	M_cmp_c = x(117);
	M_succ_c = x(118);
	M_asp_L_c = x(119);
	M_gly_L_c = x(120);
	M_ile_L_c = x(121);
	M_asn_L_c = x(122);
	M_cys_L_c = x(123);
	M_lys_L_c = x(124);
	M_his_L_c = x(125);
	M_ala_L_c = x(126);
	M_phe_L_c = x(127);
	M_pro_L_c = x(128);
	M_ser_L_c = x(129);
	M_thr_L_c = x(130);
	M_trp_L_c = x(131);
	M_tyr_L_c = x(132);
	M_val_L_c = x(133);
	M_met_L_c = x(134);
	M_leu_L_c = x(135);
	M_glu_L_c = x(136);
	M_gln_L_c = x(137);
	M_o2_c = x(138);
	M_co2_c = x(139);
	M_pi_c = x(140);
	M_nh4_c = x(141);
	M_so4_c = x(142);
	M_h_c = x(143);
	M_h2o_c = x(144);
	M_nad_c = x(145);
	M_nadh_c = x(146);
	M_nadp_c = x(147);
	M_nadph_c = x(148);
	E_R_glk_atp = x(149);
	E_R_pgi = x(150);
	E_R_pgi_R = x(151);
	E_R_pfk = x(152);
	E_R_fdp = x(153);
	E_R_fbaA = x(154);
	E_R_fbaA_R = x(155);
	E_R_tpiA = x(156);
	E_R_tpiA_R = x(157);
	E_R_gapA = x(158);
	E_R_gapA_R = x(159);
	E_R_pgk = x(160);
	E_R_pgk_R = x(161);
	E_R_gpm = x(162);
	E_R_gpm_R = x(163);
	E_R_eno = x(164);
	E_R_eno_R = x(165);
	E_R_pyk = x(166);
	E_R_pck = x(167);
	E_R_ppc = x(168);
	E_R_pdh = x(169);
	E_R_pps = x(170);
	E_R_zwf = x(171);
	E_R_zwf_R = x(172);
	E_R_pgl = x(173);
	E_R_gnd = x(174);
	E_R_rpe = x(175);
	E_R_rpe_R = x(176);
	E_R_rpi = x(177);
	E_R_rpi_R = x(178);
	E_R_talAB = x(179);
	E_R_talAB_R = x(180);
	E_R_tkt1 = x(181);
	E_R_tkt1_R = x(182);
	E_R_tkt2 = x(183);
	E_R_tkt2_R = x(184);
	E_R_edd = x(185);
	E_R_eda = x(186);
	E_R_gltA = x(187);
	E_R_acn = x(188);
	E_R_acn_R = x(189);
	E_R_icd = x(190);
	E_R_icd_R = x(191);
	E_R_sucAB = x(192);
	E_R_sucCD = x(193);
	E_R_sdh = x(194);
	E_R_frd = x(195);
	E_R_fum = x(196);
	E_R_fum_R = x(197);
	E_R_mdh = x(198);
	E_R_mdh_R = x(199);
	E_R_cyd = x(200);
	E_R_cyo = x(201);
	E_R_app = x(202);
	E_R_atp = x(203);
	E_R_nuo = x(204);
	E_R_pnt1 = x(205);
	E_R_pnt2 = x(206);
	E_R_ndh1 = x(207);
	E_R_ndh2 = x(208);
	E_R_ppa = x(209);
	E_R_aceA = x(210);
	E_R_aceB = x(211);
	E_R_maeA = x(212);
	E_R_maeB = x(213);
	E_R_pta = x(214);
	E_R_pta_R = x(215);
	E_R_ackA = x(216);
	E_R_ackA_R = x(217);
	E_R_acs = x(218);
	E_R_adhE = x(219);
	E_R_adhE_R = x(220);
	E_R_ldh = x(221);
	E_R_ldh_R = x(222);
	E_R_pflAB = x(223);
	E_R_alaAC = x(224);
	E_R_alaAC_R = x(225);
	E_R_arg = x(226);
	E_R_aspC = x(227);
	E_R_asnB = x(228);
	E_R_asnA = x(229);
	E_R_cysEMK = x(230);
	E_R_gltBD = x(231);
	E_R_gdhA = x(232);
	E_R_gdhA_R = x(233);
	E_R_glnA = x(234);
	E_R_glyA = x(235);
	E_R_his = x(236);
	E_R_ile = x(237);
	E_R_leu = x(238);
	E_R_lys = x(239);
	E_R_met = x(240);
	E_R_phe = x(241);
	E_R_pro = x(242);
	E_R_serABC = x(243);
	E_R_thr = x(244);
	E_R_trp = x(245);
	E_R_tyr = x(246);
	E_R_val = x(247);
	E_R_arg_deg = x(248);
	E_R_asp_deg = x(249);
	E_R_asn_deg = x(250);
	E_R_gly_deg = x(251);
	E_R_mglx_deg = x(252);
	E_R_ser_deg = x(253);
	E_R_pro_deg = x(254);
	E_R_thr_deg1 = x(255);
	E_R_thr_deg2 = x(256);
	E_R_thr_deg3 = x(257);
	E_R_trp_deg = x(258);
	E_R_cys_deg = x(259);
	E_R_lys_deg = x(260);
	E_R_gln_deg = x(261);
	E_R_glu_deg = x(262);
	E_R_gaba_deg1 = x(263);
	E_R_gaba_deg2 = x(264);
	E_R_chor = x(265);
	E_R_fol_e = x(266);
	E_R_fol_1 = x(267);
	E_R_fol_2a = x(268);
	E_R_fol_2b = x(269);
	E_R_fol_3 = x(270);
	E_R_fol_4 = x(271);
	E_R_gly_fol = x(272);
	E_R_gly_fol_R = x(273);
	E_R_mthfd = x(274);
	E_R_mthfd_R = x(275);
	E_R_mthfc = x(276);
	E_R_mthfc_R = x(277);
	E_R_mthfr2a = x(278);
	E_R_mthfr2b = x(279);
	E_R_prpp_syn = x(280);
	E_R_or_syn_1 = x(281);
	E_R_or_syn_2 = x(282);
	E_R_omp_syn = x(283);
	E_R_ump_syn = x(284);
	E_R_ctp_1 = x(285);
	E_R_ctp_2 = x(286);
	E_R_A_syn_1 = x(287);
	E_R_A_syn_2 = x(288);
	E_R_A_syn_3 = x(289);
	E_R_A_syn_4 = x(290);
	E_R_A_syn_5 = x(291);
	E_R_A_syn_6 = x(292);
	E_R_A_syn_7 = x(293);
	E_R_A_syn_8 = x(294);
	E_R_A_syn_9 = x(295);
	E_R_A_syn_10 = x(296);
	E_R_A_syn_12 = x(297);
	E_R_xmp_syn = x(298);
	E_R_gmp_syn = x(299);
	E_R_atp_amp = x(300);
	E_R_utp_ump = x(301);
	E_R_ctp_cmp = x(302);
	E_R_gtp_gmp = x(303);
	E_R_atp_adp = x(304);
	E_R_utp_adp = x(305);
	E_R_ctp_adp = x(306);
	E_R_gtp_adp = x(307);
	E_R_udp_utp = x(308);
	E_R_cdp_ctp = x(309);
	E_R_gdp_gtp = x(310);
	E_R_atp_ump = x(311);
	E_R_atp_cmp = x(312);
	E_R_atp_gmp = x(313);
	E_R_adk_atp = x(314);
	E_Import_o2 = x(315);
	E_Import_co2 = x(316);
	E_Import_pi = x(317);
	E_Import_nh4 = x(318);
	E_Import_so4 = x(319);
	E_Import_h2o = x(320);
	E_Export_o2 = x(321);
	E_Export_co2 = x(322);
	E_Export_pi = x(323);
	E_Export_nh4 = x(324);
	E_Export_so4 = x(325);
	E_Export_h2o = x(326);
	E_Proton_gradient = x(327);
	E_transcriptional_initiation_CAT = x(328);
	E_transcription_CAT = x(329);
	E_mRNA_degradation_CAT = x(330);
	E_translation_initiation_CAT = x(331);
	E_translation_CAT = x(332);
	E_tRNA_charging_M_ala_L_c_CAT = x(333);
	E_tRNA_charging_M_arg_L_c_CAT = x(334);
	E_tRNA_charging_M_asn_L_c_CAT = x(335);
	E_tRNA_charging_M_asp_L_c_CAT = x(336);
	E_tRNA_charging_M_cys_L_c_CAT = x(337);
	E_tRNA_charging_M_glu_L_c_CAT = x(338);
	E_tRNA_charging_M_gln_L_c_CAT = x(339);
	E_tRNA_charging_M_gly_L_c_CAT = x(340);
	E_tRNA_charging_M_his_L_c_CAT = x(341);
	E_tRNA_charging_M_ile_L_c_CAT = x(342);
	E_tRNA_charging_M_leu_L_c_CAT = x(343);
	E_tRNA_charging_M_lys_L_c_CAT = x(344);
	E_tRNA_charging_M_met_L_c_CAT = x(345);
	E_tRNA_charging_M_phe_L_c_CAT = x(346);
	E_tRNA_charging_M_pro_L_c_CAT = x(347);
	E_tRNA_charging_M_ser_L_c_CAT = x(348);
	E_tRNA_charging_M_thr_L_c_CAT = x(349);
	E_tRNA_charging_M_trp_L_c_CAT = x(350);
	E_tRNA_charging_M_tyr_L_c_CAT = x(351);
	E_tRNA_charging_M_val_L_c_CAT = x(352);
    
    rate_constant_array = model.K;
    saturation_constant_array = model.saturation_constants;
    
    % Formulate the kinetic rate vector - 
	%rate_constant_array = data_dictionary("RATE_CONSTANT_ARRAY");
	%saturation_constant_array = data_dictionary("SATURATION_CONSTANT_ARRAY");
	rate_vector = [];
	

	% R_glk_atp: M_atp_c+M_glc_D_c = M_adp_c+M_g6p_c+M_h_c
	tmp_reaction = rate_constant_array(1)*(E_R_glk_atp)*(M_atp_c/(saturation_constant_array(1,106) + M_atp_c))*(M_glc_D_c/(saturation_constant_array(1,100) + M_glc_D_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pgi: M_g6p_c = M_f6p_c
	tmp_reaction = rate_constant_array(2)*(E_R_pgi)*(M_g6p_c/(saturation_constant_array(2,1) + M_g6p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pgi_R: M_f6p_c = M_g6p_c
	tmp_reaction = rate_constant_array(3)*(E_R_pgi_R)*(M_f6p_c/(saturation_constant_array(3,2) + M_f6p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pfk: M_atp_c+M_f6p_c = M_adp_c+M_fdp_c
	tmp_reaction = rate_constant_array(4)*(E_R_pfk)*(M_atp_c/(saturation_constant_array(4,106) + M_atp_c))*(M_f6p_c/(saturation_constant_array(4,2) + M_f6p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fdp: M_fdp_c+M_h2o_c = M_f6p_c+M_pi_c
	tmp_reaction = rate_constant_array(5)*(E_R_fdp)*(M_fdp_c/(saturation_constant_array(5,3) + M_fdp_c))*(M_h2o_c/(saturation_constant_array(5,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fbaA: M_fdp_c = M_dhap_c+M_g3p_c
	tmp_reaction = rate_constant_array(6)*(E_R_fbaA)*(M_fdp_c/(saturation_constant_array(6,3) + M_fdp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fbaA_R: M_dhap_c+M_g3p_c = M_fdp_c
	tmp_reaction = rate_constant_array(7)*(E_R_fbaA_R)*(M_dhap_c/(saturation_constant_array(7,4) + M_dhap_c))*(M_g3p_c/(saturation_constant_array(7,5) + M_g3p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_tpiA: M_dhap_c = M_g3p_c
	tmp_reaction = rate_constant_array(8)*(E_R_tpiA)*(M_dhap_c/(saturation_constant_array(8,4) + M_dhap_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_tpiA_R: M_g3p_c = M_dhap_c
	tmp_reaction = rate_constant_array(9)*(E_R_tpiA_R)*(M_g3p_c/(saturation_constant_array(9,5) + M_g3p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gapA: M_g3p_c+M_nad_c+M_pi_c = M_13dpg_c+M_h_c+M_nadh_c
	tmp_reaction = rate_constant_array(10)*(E_R_gapA)*(M_g3p_c/(saturation_constant_array(10,5) + M_g3p_c))*(M_nad_c/(saturation_constant_array(10,145) + M_nad_c))*(M_pi_c/(saturation_constant_array(10,140) + M_pi_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gapA_R: M_13dpg_c+M_h_c+M_nadh_c = M_g3p_c+M_nad_c+M_pi_c
	tmp_reaction = rate_constant_array(11)*(E_R_gapA_R)*(M_13dpg_c/(saturation_constant_array(11,6) + M_13dpg_c))*(M_h_c/(saturation_constant_array(11,143) + M_h_c))*(M_nadh_c/(saturation_constant_array(11,146) + M_nadh_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pgk: M_13dpg_c+M_adp_c = M_3pg_c+M_atp_c
	tmp_reaction = rate_constant_array(12)*(E_R_pgk)*(M_13dpg_c/(saturation_constant_array(12,6) + M_13dpg_c))*(M_adp_c/(saturation_constant_array(12,107) + M_adp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pgk_R: M_3pg_c+M_atp_c = M_13dpg_c+M_adp_c
	tmp_reaction = rate_constant_array(13)*(E_R_pgk_R)*(M_3pg_c/(saturation_constant_array(13,7) + M_3pg_c))*(M_atp_c/(saturation_constant_array(13,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gpm: M_3pg_c = M_2pg_c
	tmp_reaction = rate_constant_array(14)*(E_R_gpm)*(M_3pg_c/(saturation_constant_array(14,7) + M_3pg_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gpm_R: M_2pg_c = M_3pg_c
	tmp_reaction = rate_constant_array(15)*(E_R_gpm_R)*(M_2pg_c/(saturation_constant_array(15,8) + M_2pg_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_eno: M_2pg_c = M_h2o_c+M_pep_c
	tmp_reaction = rate_constant_array(16)*(E_R_eno)*(M_2pg_c/(saturation_constant_array(16,8) + M_2pg_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_eno_R: M_h2o_c+M_pep_c = M_2pg_c
	tmp_reaction = rate_constant_array(17)*(E_R_eno_R)*(M_h2o_c/(saturation_constant_array(17,144) + M_h2o_c))*(M_pep_c/(saturation_constant_array(17,101) + M_pep_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pyk: M_adp_c+M_pep_c = M_atp_c+M_pyr_c
	tmp_reaction = rate_constant_array(18)*(E_R_pyk)*(M_adp_c/(saturation_constant_array(18,107) + M_adp_c))*(M_pep_c/(saturation_constant_array(18,101) + M_pep_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pck: M_atp_c+M_oaa_c = M_adp_c+M_co2_c+M_pep_c
	tmp_reaction = rate_constant_array(19)*(E_R_pck)*(M_atp_c/(saturation_constant_array(19,106) + M_atp_c))*(M_oaa_c/(saturation_constant_array(19,9) + M_oaa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ppc: M_co2_c+M_h2o_c+M_pep_c = M_oaa_c+M_pi_c
	tmp_reaction = rate_constant_array(20)*(E_R_ppc)*(M_co2_c/(saturation_constant_array(20,139) + M_co2_c))*(M_h2o_c/(saturation_constant_array(20,144) + M_h2o_c))*(M_pep_c/(saturation_constant_array(20,101) + M_pep_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pdh: M_coa_c+M_nad_c+M_pyr_c = M_accoa_c+M_co2_c+M_nadh_c+M_h_c
	tmp_reaction = rate_constant_array(21)*(E_R_pdh)*(M_coa_c/(saturation_constant_array(21,10) + M_coa_c))*(M_nad_c/(saturation_constant_array(21,145) + M_nad_c))*(M_pyr_c/(saturation_constant_array(21,102) + M_pyr_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pps: M_atp_c+M_h2o_c+M_pyr_c = M_amp_c+M_pep_c+M_pi_c
	tmp_reaction = rate_constant_array(22)*(E_R_pps)*(M_atp_c/(saturation_constant_array(22,106) + M_atp_c))*(M_h2o_c/(saturation_constant_array(22,144) + M_h2o_c))*(M_pyr_c/(saturation_constant_array(22,102) + M_pyr_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_zwf: M_g6p_c+M_nadp_c = M_6pgl_c+M_h_c+M_nadph_c
	tmp_reaction = rate_constant_array(23)*(E_R_zwf)*(M_g6p_c/(saturation_constant_array(23,1) + M_g6p_c))*(M_nadp_c/(saturation_constant_array(23,147) + M_nadp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_zwf_R: M_6pgl_c+M_h_c+M_nadph_c = M_g6p_c+M_nadp_c
	tmp_reaction = rate_constant_array(24)*(E_R_zwf_R)*(M_6pgl_c/(saturation_constant_array(24,12) + M_6pgl_c))*(M_h_c/(saturation_constant_array(24,143) + M_h_c))*(M_nadph_c/(saturation_constant_array(24,148) + M_nadph_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pgl: M_6pgl_c+M_h2o_c = M_6pgc_c
	tmp_reaction = rate_constant_array(25)*(E_R_pgl)*(M_6pgl_c/(saturation_constant_array(25,12) + M_6pgl_c))*(M_h2o_c/(saturation_constant_array(25,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gnd: M_6pgc_c+M_nadp_c = M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c
	tmp_reaction = rate_constant_array(26)*(E_R_gnd)*(M_6pgc_c/(saturation_constant_array(26,13) + M_6pgc_c))*(M_nadp_c/(saturation_constant_array(26,147) + M_nadp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_rpe: M_ru5p_D_c = M_xu5p_D_c
	tmp_reaction = rate_constant_array(27)*(E_R_rpe)*(M_ru5p_D_c/(saturation_constant_array(27,14) + M_ru5p_D_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_rpe_R: M_xu5p_D_c = M_ru5p_D_c
	tmp_reaction = rate_constant_array(28)*(E_R_rpe_R)*(M_xu5p_D_c/(saturation_constant_array(28,15) + M_xu5p_D_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_rpi: M_r5p_c = M_ru5p_D_c
	tmp_reaction = rate_constant_array(29)*(E_R_rpi)*(M_r5p_c/(saturation_constant_array(29,16) + M_r5p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_rpi_R: M_ru5p_D_c = M_r5p_c
	tmp_reaction = rate_constant_array(30)*(E_R_rpi_R)*(M_ru5p_D_c/(saturation_constant_array(30,14) + M_ru5p_D_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_talAB: M_g3p_c+M_s7p_c = M_e4p_c+M_f6p_c
	tmp_reaction = rate_constant_array(31)*(E_R_talAB)*(M_g3p_c/(saturation_constant_array(31,5) + M_g3p_c))*(M_s7p_c/(saturation_constant_array(31,17) + M_s7p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_talAB_R: M_e4p_c+M_f6p_c = M_g3p_c+M_s7p_c
	tmp_reaction = rate_constant_array(32)*(E_R_talAB_R)*(M_e4p_c/(saturation_constant_array(32,18) + M_e4p_c))*(M_f6p_c/(saturation_constant_array(32,2) + M_f6p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_tkt1: M_r5p_c+M_xu5p_D_c = M_g3p_c+M_s7p_c
	tmp_reaction = rate_constant_array(33)*(E_R_tkt1)*(M_r5p_c/(saturation_constant_array(33,16) + M_r5p_c))*(M_xu5p_D_c/(saturation_constant_array(33,15) + M_xu5p_D_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_tkt1_R: M_g3p_c+M_s7p_c = M_r5p_c+M_xu5p_D_c
	tmp_reaction = rate_constant_array(34)*(E_R_tkt1_R)*(M_g3p_c/(saturation_constant_array(34,5) + M_g3p_c))*(M_s7p_c/(saturation_constant_array(34,17) + M_s7p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_tkt2: M_e4p_c+M_xu5p_D_c = M_f6p_c+M_g3p_c
	tmp_reaction = rate_constant_array(35)*(E_R_tkt2)*(M_e4p_c/(saturation_constant_array(35,18) + M_e4p_c))*(M_xu5p_D_c/(saturation_constant_array(35,15) + M_xu5p_D_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_tkt2_R: M_f6p_c+M_g3p_c = M_e4p_c+M_xu5p_D_c
	tmp_reaction = rate_constant_array(36)*(E_R_tkt2_R)*(M_f6p_c/(saturation_constant_array(36,2) + M_f6p_c))*(M_g3p_c/(saturation_constant_array(36,5) + M_g3p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_edd: M_6pgc_c = M_2ddg6p_c+M_h2o_c
	tmp_reaction = rate_constant_array(37)*(E_R_edd)*(M_6pgc_c/(saturation_constant_array(37,13) + M_6pgc_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_eda: M_2ddg6p_c = M_g3p_c+M_pyr_c
	tmp_reaction = rate_constant_array(38)*(E_R_eda)*(M_2ddg6p_c/(saturation_constant_array(38,19) + M_2ddg6p_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gltA: M_accoa_c+M_h2o_c+M_oaa_c = M_cit_c+M_coa_c
	tmp_reaction = rate_constant_array(39)*(E_R_gltA)*(M_accoa_c/(saturation_constant_array(39,11) + M_accoa_c))*(M_h2o_c/(saturation_constant_array(39,144) + M_h2o_c))*(M_oaa_c/(saturation_constant_array(39,9) + M_oaa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_acn: M_cit_c = M_icit_c
	tmp_reaction = rate_constant_array(40)*(E_R_acn)*(M_cit_c/(saturation_constant_array(40,20) + M_cit_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_acn_R: M_icit_c = M_cit_c
	tmp_reaction = rate_constant_array(41)*(E_R_acn_R)*(M_icit_c/(saturation_constant_array(41,21) + M_icit_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_icd: M_icit_c+M_nadp_c = M_akg_c+M_co2_c+M_nadph_c+M_h_c
	tmp_reaction = rate_constant_array(42)*(E_R_icd)*(M_icit_c/(saturation_constant_array(42,21) + M_icit_c))*(M_nadp_c/(saturation_constant_array(42,147) + M_nadp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_icd_R: M_akg_c+M_co2_c+M_nadph_c+M_h_c = M_icit_c+M_nadp_c
	tmp_reaction = rate_constant_array(43)*(E_R_icd_R)*(M_akg_c/(saturation_constant_array(43,22) + M_akg_c))*(M_co2_c/(saturation_constant_array(43,139) + M_co2_c))*(M_nadph_c/(saturation_constant_array(43,148) + M_nadph_c))*(M_h_c/(saturation_constant_array(43,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_sucAB: M_akg_c+M_coa_c+M_nad_c = M_co2_c+M_nadh_c+M_succoa_c+M_h_c
	tmp_reaction = rate_constant_array(44)*(E_R_sucAB)*(M_akg_c/(saturation_constant_array(44,22) + M_akg_c))*(M_coa_c/(saturation_constant_array(44,10) + M_coa_c))*(M_nad_c/(saturation_constant_array(44,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_sucCD: M_adp_c+M_pi_c+M_succoa_c = M_atp_c+M_coa_c+M_succ_c
	tmp_reaction = rate_constant_array(45)*(E_R_sucCD)*(M_adp_c/(saturation_constant_array(45,107) + M_adp_c))*(M_pi_c/(saturation_constant_array(45,140) + M_pi_c))*(M_succoa_c/(saturation_constant_array(45,23) + M_succoa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_sdh: M_q8_c+M_succ_c = M_fum_c+M_q8h2_c
	tmp_reaction = rate_constant_array(46)*(E_R_sdh)*(M_q8_c/(saturation_constant_array(46,24) + M_q8_c))*(M_succ_c/(saturation_constant_array(46,118) + M_succ_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_frd: M_fum_c+M_mql8_c = M_mqn8_c+M_succ_c
	tmp_reaction = rate_constant_array(47)*(E_R_frd)*(M_fum_c/(saturation_constant_array(47,25) + M_fum_c))*(M_mql8_c/(saturation_constant_array(47,27) + M_mql8_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fum: M_fum_c+M_h2o_c = M_mal_L_c
	tmp_reaction = rate_constant_array(48)*(E_R_fum)*(M_fum_c/(saturation_constant_array(48,25) + M_fum_c))*(M_h2o_c/(saturation_constant_array(48,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fum_R: M_mal_L_c = M_fum_c+M_h2o_c
	tmp_reaction = rate_constant_array(49)*(E_R_fum_R)*(M_mal_L_c/(saturation_constant_array(49,105) + M_mal_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mdh: M_mal_L_c+M_nad_c = M_oaa_c+M_h_c+M_nadh_c
	tmp_reaction = rate_constant_array(50)*(E_R_mdh)*(M_mal_L_c/(saturation_constant_array(50,105) + M_mal_L_c))*(M_nad_c/(saturation_constant_array(50,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mdh_R: M_oaa_c+M_h_c+M_nadh_c = M_mal_L_c+M_nad_c
	tmp_reaction = rate_constant_array(51)*(E_R_mdh_R)*(M_oaa_c/(saturation_constant_array(51,9) + M_oaa_c))*(M_h_c/(saturation_constant_array(51,143) + M_h_c))*(M_nadh_c/(saturation_constant_array(51,146) + M_nadh_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_cyd: 2*M_h_c+0.5*M_o2_c+M_q8h2_c = M_h2o_c+M_q8_c+2*M_h_e
	tmp_reaction = rate_constant_array(52)*(E_R_cyd)*(M_h_c/(saturation_constant_array(52,143) + M_h_c))*(M_o2_c/(saturation_constant_array(52,138) + M_o2_c))*(M_q8h2_c/(saturation_constant_array(52,26) + M_q8h2_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_cyo: 4*M_h_c+0.5*M_o2_c+M_q8h2_c = M_h2o_c+M_q8_c+4*M_h_e
	tmp_reaction = rate_constant_array(53)*(E_R_cyo)*(M_h_c/(saturation_constant_array(53,143) + M_h_c))*(M_o2_c/(saturation_constant_array(53,138) + M_o2_c))*(M_q8h2_c/(saturation_constant_array(53,26) + M_q8h2_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_app: 2*M_h_c+M_mql8_c+0.5*M_o2_c = M_h2o_c+M_mqn8_c+2*M_h_e
	tmp_reaction = rate_constant_array(54)*(E_R_app)*(M_h_c/(saturation_constant_array(54,143) + M_h_c))*(M_mql8_c/(saturation_constant_array(54,27) + M_mql8_c))*(M_o2_c/(saturation_constant_array(54,138) + M_o2_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_atp: M_adp_c+M_pi_c+4*M_h_e = M_atp_c+4*M_h_c+M_h2o_c
	tmp_reaction = rate_constant_array(55)*(E_R_atp)*(M_adp_c/(saturation_constant_array(55,107) + M_adp_c))*(M_pi_c/(saturation_constant_array(55,140) + M_pi_c))*(M_h_e/(saturation_constant_array(55,29) + M_h_e));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_nuo: 3*M_h_c+M_nadh_c+M_q8_c = M_nad_c+M_q8h2_c+2*M_h_e
	tmp_reaction = rate_constant_array(56)*(E_R_nuo)*(M_h_c/(saturation_constant_array(56,143) + M_h_c))*(M_nadh_c/(saturation_constant_array(56,146) + M_nadh_c))*(M_q8_c/(saturation_constant_array(56,24) + M_q8_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pnt1: M_nad_c+M_nadph_c = M_nadh_c+M_nadp_c
	tmp_reaction = rate_constant_array(57)*(E_R_pnt1)*(M_nad_c/(saturation_constant_array(57,145) + M_nad_c))*(M_nadph_c/(saturation_constant_array(57,148) + M_nadph_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pnt2: M_nadh_c+M_nadp_c+2*M_h_e = 2*M_h_c+M_nad_c+M_nadph_c
	tmp_reaction = rate_constant_array(58)*(E_R_pnt2)*(M_nadh_c/(saturation_constant_array(58,146) + M_nadh_c))*(M_nadp_c/(saturation_constant_array(58,147) + M_nadp_c))*(M_h_e/(saturation_constant_array(58,29) + M_h_e));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ndh1: M_h_c+M_nadh_c+M_q8_c = M_nad_c+M_q8h2_c
	tmp_reaction = rate_constant_array(59)*(E_R_ndh1)*(M_h_c/(saturation_constant_array(59,143) + M_h_c))*(M_nadh_c/(saturation_constant_array(59,146) + M_nadh_c))*(M_q8_c/(saturation_constant_array(59,24) + M_q8_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ndh2: M_h_c+M_mqn8_c+M_nadh_c = M_mql8_c+M_nad_c
	tmp_reaction = rate_constant_array(60)*(E_R_ndh2)*(M_h_c/(saturation_constant_array(60,143) + M_h_c))*(M_mqn8_c/(saturation_constant_array(60,28) + M_mqn8_c))*(M_nadh_c/(saturation_constant_array(60,146) + M_nadh_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ppa: M_ppi_c+M_h2o_c = 2*M_pi_c
	tmp_reaction = rate_constant_array(61)*(E_R_ppa)*(M_ppi_c/(saturation_constant_array(61,30) + M_ppi_c))*(M_h2o_c/(saturation_constant_array(61,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_aceA: M_icit_c = M_glx_c+M_succ_c
	tmp_reaction = rate_constant_array(62)*(E_R_aceA)*(M_icit_c/(saturation_constant_array(62,21) + M_icit_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_aceB: M_accoa_c+M_glx_c+M_h2o_c = M_coa_c+M_mal_L_c
	tmp_reaction = rate_constant_array(63)*(E_R_aceB)*(M_accoa_c/(saturation_constant_array(63,11) + M_accoa_c))*(M_glx_c/(saturation_constant_array(63,31) + M_glx_c))*(M_h2o_c/(saturation_constant_array(63,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_maeA: M_mal_L_c+M_nad_c = M_co2_c+M_nadh_c+M_pyr_c+M_h_c
	tmp_reaction = rate_constant_array(64)*(E_R_maeA)*(M_mal_L_c/(saturation_constant_array(64,105) + M_mal_L_c))*(M_nad_c/(saturation_constant_array(64,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_maeB: M_mal_L_c+M_nadp_c = M_co2_c+M_nadph_c+M_pyr_c+M_h_c
	tmp_reaction = rate_constant_array(65)*(E_R_maeB)*(M_mal_L_c/(saturation_constant_array(65,105) + M_mal_L_c))*(M_nadp_c/(saturation_constant_array(65,147) + M_nadp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pta: M_accoa_c+M_pi_c = M_actp_c+M_coa_c
	tmp_reaction = rate_constant_array(66)*(E_R_pta)*(M_accoa_c/(saturation_constant_array(66,11) + M_accoa_c))*(M_pi_c/(saturation_constant_array(66,140) + M_pi_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pta_R: M_actp_c+M_coa_c = M_accoa_c+M_pi_c
	tmp_reaction = rate_constant_array(67)*(E_R_pta_R)*(M_actp_c/(saturation_constant_array(67,32) + M_actp_c))*(M_coa_c/(saturation_constant_array(67,10) + M_coa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ackA: M_actp_c+M_adp_c = M_ac_c+M_atp_c
	tmp_reaction = rate_constant_array(68)*(E_R_ackA)*(M_actp_c/(saturation_constant_array(68,32) + M_actp_c))*(M_adp_c/(saturation_constant_array(68,107) + M_adp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ackA_R: M_ac_c+M_atp_c = M_actp_c+M_adp_c
	tmp_reaction = rate_constant_array(69)*(E_R_ackA_R)*(M_ac_c/(saturation_constant_array(69,103) + M_ac_c))*(M_atp_c/(saturation_constant_array(69,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_acs: M_ac_c+M_atp_c+M_coa_c = M_accoa_c+M_amp_c+M_ppi_c
	tmp_reaction = rate_constant_array(70)*(E_R_acs)*(M_ac_c/(saturation_constant_array(70,103) + M_ac_c))*(M_atp_c/(saturation_constant_array(70,106) + M_atp_c))*(M_coa_c/(saturation_constant_array(70,10) + M_coa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_adhE: M_accoa_c+2*M_h_c+2*M_nadh_c = M_coa_c+M_etoh_c+2*M_nad_c
	tmp_reaction = rate_constant_array(71)*(E_R_adhE)*(M_accoa_c/(saturation_constant_array(71,11) + M_accoa_c))*(M_h_c/(saturation_constant_array(71,143) + M_h_c))*(M_nadh_c/(saturation_constant_array(71,146) + M_nadh_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_adhE_R: M_coa_c+M_etoh_c+2*M_nad_c = M_accoa_c+2*M_h_c+2*M_nadh_c
	tmp_reaction = rate_constant_array(72)*(E_R_adhE_R)*(M_coa_c/(saturation_constant_array(72,10) + M_coa_c))*(M_etoh_c/(saturation_constant_array(72,33) + M_etoh_c))*(M_nad_c/(saturation_constant_array(72,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ldh: M_pyr_c+M_nadh_c+M_h_c = M_lac_D_c+M_nad_c
	tmp_reaction = rate_constant_array(73)*(E_R_ldh)*(M_pyr_c/(saturation_constant_array(73,102) + M_pyr_c))*(M_nadh_c/(saturation_constant_array(73,146) + M_nadh_c))*(M_h_c/(saturation_constant_array(73,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ldh_R: M_lac_D_c+M_nad_c = M_pyr_c+M_nadh_c+M_h_c
	tmp_reaction = rate_constant_array(74)*(E_R_ldh_R)*(M_lac_D_c/(saturation_constant_array(74,104) + M_lac_D_c))*(M_nad_c/(saturation_constant_array(74,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pflAB: M_coa_c+M_pyr_c = M_accoa_c+M_for_c
	tmp_reaction = rate_constant_array(75)*(E_R_pflAB)*(M_coa_c/(saturation_constant_array(75,10) + M_coa_c))*(M_pyr_c/(saturation_constant_array(75,102) + M_pyr_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_alaAC: M_pyr_c+M_glu_L_c = M_ala_L_c+M_akg_c
	tmp_reaction = rate_constant_array(76)*(E_R_alaAC)*(M_pyr_c/(saturation_constant_array(76,102) + M_pyr_c))*(M_glu_L_c/(saturation_constant_array(76,136) + M_glu_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_alaAC_R: M_ala_L_c+M_akg_c = M_pyr_c+M_glu_L_c
	tmp_reaction = rate_constant_array(77)*(E_R_alaAC_R)*(M_ala_L_c/(saturation_constant_array(77,126) + M_ala_L_c))*(M_akg_c/(saturation_constant_array(77,22) + M_akg_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_arg: M_accoa_c+2*M_glu_L_c+3*M_atp_c+M_nadph_c+M_h_c+M_h2o_c+M_nh3_c+M_co2_c+M_asp_L_c = M_coa_c+2*M_adp_c+2*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c+M_arg_L_c
	tmp_reaction = rate_constant_array(78)*(E_R_arg)*(M_accoa_c/(saturation_constant_array(78,11) + M_accoa_c))*(M_glu_L_c/(saturation_constant_array(78,136) + M_glu_L_c))*(M_atp_c/(saturation_constant_array(78,106) + M_atp_c))*(M_nadph_c/(saturation_constant_array(78,148) + M_nadph_c))*(M_h_c/(saturation_constant_array(78,143) + M_h_c))*(M_h2o_c/(saturation_constant_array(78,144) + M_h2o_c))*(M_nh3_c/(saturation_constant_array(78,35) + M_nh3_c))*(M_co2_c/(saturation_constant_array(78,139) + M_co2_c))*(M_asp_L_c/(saturation_constant_array(78,119) + M_asp_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_aspC: M_glu_L_c+M_oaa_c = M_asp_L_c+M_akg_c
	tmp_reaction = rate_constant_array(79)*(E_R_aspC)*(M_glu_L_c/(saturation_constant_array(79,136) + M_glu_L_c))*(M_oaa_c/(saturation_constant_array(79,9) + M_oaa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_asnB: M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c = M_asn_L_c+M_glu_L_c+M_ppi_c+M_amp_c
	tmp_reaction = rate_constant_array(80)*(E_R_asnB)*(M_asp_L_c/(saturation_constant_array(80,119) + M_asp_L_c))*(M_gln_L_c/(saturation_constant_array(80,137) + M_gln_L_c))*(M_h2o_c/(saturation_constant_array(80,144) + M_h2o_c))*(M_atp_c/(saturation_constant_array(80,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_asnA: M_asp_L_c+M_atp_c+M_nh3_c = M_asn_L_c+M_ppi_c+M_amp_c
	tmp_reaction = rate_constant_array(81)*(E_R_asnA)*(M_asp_L_c/(saturation_constant_array(81,119) + M_asp_L_c))*(M_atp_c/(saturation_constant_array(81,106) + M_atp_c))*(M_nh3_c/(saturation_constant_array(81,35) + M_nh3_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_cysEMK: M_ser_L_c+M_accoa_c+M_h2s_c = M_cys_L_c+M_coa_c+M_ac_c
	tmp_reaction = rate_constant_array(82)*(E_R_cysEMK)*(M_ser_L_c/(saturation_constant_array(82,129) + M_ser_L_c))*(M_accoa_c/(saturation_constant_array(82,11) + M_accoa_c))*(M_h2s_c/(saturation_constant_array(82,37) + M_h2s_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gltBD: M_gln_L_c+M_akg_c+M_nadph_c+M_h_c = 2.0*M_glu_L_c+M_nadp_c
	tmp_reaction = rate_constant_array(83)*(E_R_gltBD)*(M_gln_L_c/(saturation_constant_array(83,137) + M_gln_L_c))*(M_akg_c/(saturation_constant_array(83,22) + M_akg_c))*(M_nadph_c/(saturation_constant_array(83,148) + M_nadph_c))*(M_h_c/(saturation_constant_array(83,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gdhA: M_akg_c+M_nadph_c+M_nh3_c+M_h_c = M_glu_L_c+M_h2o_c+M_nadp_c
	tmp_reaction = rate_constant_array(84)*(E_R_gdhA)*(M_akg_c/(saturation_constant_array(84,22) + M_akg_c))*(M_nadph_c/(saturation_constant_array(84,148) + M_nadph_c))*(M_nh3_c/(saturation_constant_array(84,35) + M_nh3_c))*(M_h_c/(saturation_constant_array(84,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gdhA_R: M_glu_L_c+M_h2o_c+M_nadp_c = M_akg_c+M_nadph_c+M_nh3_c+M_h_c
	tmp_reaction = rate_constant_array(85)*(E_R_gdhA_R)*(M_glu_L_c/(saturation_constant_array(85,136) + M_glu_L_c))*(M_h2o_c/(saturation_constant_array(85,144) + M_h2o_c))*(M_nadp_c/(saturation_constant_array(85,147) + M_nadp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_glnA: M_glu_L_c+M_atp_c+M_nh3_c = M_gln_L_c+M_adp_c+M_pi_c
	tmp_reaction = rate_constant_array(86)*(E_R_glnA)*(M_glu_L_c/(saturation_constant_array(86,136) + M_glu_L_c))*(M_atp_c/(saturation_constant_array(86,106) + M_atp_c))*(M_nh3_c/(saturation_constant_array(86,35) + M_nh3_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_glyA: M_ser_L_c+M_thf_c = M_gly_L_c+M_h2o_c+M_mlthf_c
	tmp_reaction = rate_constant_array(87)*(E_R_glyA)*(M_ser_L_c/(saturation_constant_array(87,129) + M_ser_L_c))*(M_thf_c/(saturation_constant_array(87,38) + M_thf_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_his: M_gln_L_c+M_r5p_c+2.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c = M_his_L_c+M_akg_c+M_aicar_c+2.0*M_nadh_c+M_amp_c+M_pi_c+2.0*M_ppi_c+2.0*M_h_c
	tmp_reaction = rate_constant_array(88)*(E_R_his)*(M_gln_L_c/(saturation_constant_array(88,137) + M_gln_L_c))*(M_r5p_c/(saturation_constant_array(88,16) + M_r5p_c))*(M_atp_c/(saturation_constant_array(88,106) + M_atp_c))*(M_nad_c/(saturation_constant_array(88,145) + M_nad_c))*(M_h2o_c/(saturation_constant_array(88,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ile: M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c = M_ile_L_c+M_h2o_c+M_nh3_c+M_co2_c+M_nadp_c+M_akg_c
	tmp_reaction = rate_constant_array(89)*(E_R_ile)*(M_thr_L_c/(saturation_constant_array(89,130) + M_thr_L_c))*(M_h_c/(saturation_constant_array(89,143) + M_h_c))*(M_pyr_c/(saturation_constant_array(89,102) + M_pyr_c))*(M_nadph_c/(saturation_constant_array(89,148) + M_nadph_c))*(M_glu_L_c/(saturation_constant_array(89,136) + M_glu_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_leu: 2.0*M_pyr_c+M_glu_L_c+M_nad_c+M_nadph_c+M_accoa_c = M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
	tmp_reaction = rate_constant_array(90)*(E_R_leu)*(M_pyr_c/(saturation_constant_array(90,102) + M_pyr_c))*(M_glu_L_c/(saturation_constant_array(90,136) + M_glu_L_c))*(M_nad_c/(saturation_constant_array(90,145) + M_nad_c))*(M_nadph_c/(saturation_constant_array(90,148) + M_nadph_c))*(M_accoa_c/(saturation_constant_array(90,11) + M_accoa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_lys: M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c = M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
	tmp_reaction = rate_constant_array(91)*(E_R_lys)*(M_asp_L_c/(saturation_constant_array(91,119) + M_asp_L_c))*(M_atp_c/(saturation_constant_array(91,106) + M_atp_c))*(M_nadph_c/(saturation_constant_array(91,148) + M_nadph_c))*(M_h_c/(saturation_constant_array(91,143) + M_h_c))*(M_pyr_c/(saturation_constant_array(91,102) + M_pyr_c))*(M_succoa_c/(saturation_constant_array(91,23) + M_succoa_c))*(M_glu_L_c/(saturation_constant_array(91,136) + M_glu_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_met: M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c+M_h2o_c+2.0*M_h_c = M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh3_c+M_pyr_c
	tmp_reaction = rate_constant_array(92)*(E_R_met)*(M_asp_L_c/(saturation_constant_array(92,119) + M_asp_L_c))*(M_cys_L_c/(saturation_constant_array(92,123) + M_cys_L_c))*(M_succoa_c/(saturation_constant_array(92,23) + M_succoa_c))*(M_atp_c/(saturation_constant_array(92,106) + M_atp_c))*(M_nadph_c/(saturation_constant_array(92,148) + M_nadph_c))*(M_5mthf_c/(saturation_constant_array(92,41) + M_5mthf_c))*(M_h2o_c/(saturation_constant_array(92,144) + M_h2o_c))*(M_h_c/(saturation_constant_array(92,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_phe: M_chor_c+M_glu_L_c = M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
	tmp_reaction = rate_constant_array(93)*(E_R_phe)*(M_chor_c/(saturation_constant_array(93,42) + M_chor_c))*(M_glu_L_c/(saturation_constant_array(93,136) + M_glu_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pro: M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c = M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c
	tmp_reaction = rate_constant_array(94)*(E_R_pro)*(M_glu_L_c/(saturation_constant_array(94,136) + M_glu_L_c))*(M_atp_c/(saturation_constant_array(94,106) + M_atp_c))*(M_h_c/(saturation_constant_array(94,143) + M_h_c))*(M_nadph_c/(saturation_constant_array(94,148) + M_nadph_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_serABC: M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c = M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
	tmp_reaction = rate_constant_array(95)*(E_R_serABC)*(M_3pg_c/(saturation_constant_array(95,7) + M_3pg_c))*(M_nad_c/(saturation_constant_array(95,145) + M_nad_c))*(M_glu_L_c/(saturation_constant_array(95,136) + M_glu_L_c))*(M_h2o_c/(saturation_constant_array(95,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_thr: M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_h2o_c = M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c
	tmp_reaction = rate_constant_array(96)*(E_R_thr)*(M_asp_L_c/(saturation_constant_array(96,119) + M_asp_L_c))*(M_atp_c/(saturation_constant_array(96,106) + M_atp_c))*(M_nadph_c/(saturation_constant_array(96,148) + M_nadph_c))*(M_h_c/(saturation_constant_array(96,143) + M_h_c))*(M_h2o_c/(saturation_constant_array(96,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_trp: M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+M_atp_c = M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+M_amp_c
	tmp_reaction = rate_constant_array(97)*(E_R_trp)*(M_chor_c/(saturation_constant_array(97,42) + M_chor_c))*(M_gln_L_c/(saturation_constant_array(97,137) + M_gln_L_c))*(M_ser_L_c/(saturation_constant_array(97,129) + M_ser_L_c))*(M_r5p_c/(saturation_constant_array(97,16) + M_r5p_c))*(M_atp_c/(saturation_constant_array(97,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_tyr: M_chor_c+M_glu_L_c+M_nad_c = M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c+M_h_c
	tmp_reaction = rate_constant_array(98)*(E_R_tyr)*(M_chor_c/(saturation_constant_array(98,42) + M_chor_c))*(M_glu_L_c/(saturation_constant_array(98,136) + M_glu_L_c))*(M_nad_c/(saturation_constant_array(98,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_val: 2.0*M_pyr_c+M_h_c+M_nadph_c+M_glu_L_c = M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
	tmp_reaction = rate_constant_array(99)*(E_R_val)*(M_pyr_c/(saturation_constant_array(99,102) + M_pyr_c))*(M_h_c/(saturation_constant_array(99,143) + M_h_c))*(M_nadph_c/(saturation_constant_array(99,148) + M_nadph_c))*(M_glu_L_c/(saturation_constant_array(99,136) + M_glu_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_arg_deg: M_arg_L_c+4.0*M_h2o_c+M_nad_c+M_akg_c+M_succoa_c = M_h_c+M_co2_c+2.0*M_glu_L_c+2.0*M_nh3_c+M_nadh_c+M_succ_c+M_coa_c
	tmp_reaction = rate_constant_array(100)*(E_R_arg_deg)*(M_arg_L_c/(saturation_constant_array(100,36) + M_arg_L_c))*(M_h2o_c/(saturation_constant_array(100,144) + M_h2o_c))*(M_nad_c/(saturation_constant_array(100,145) + M_nad_c))*(M_akg_c/(saturation_constant_array(100,22) + M_akg_c))*(M_succoa_c/(saturation_constant_array(100,23) + M_succoa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_asp_deg: M_asp_L_c = M_fum_c+M_nh3_c
	tmp_reaction = rate_constant_array(101)*(E_R_asp_deg)*(M_asp_L_c/(saturation_constant_array(101,119) + M_asp_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_asn_deg: M_asn_L_c+M_amp_c+M_ppi_c = M_nh3_c+M_asp_L_c+M_atp_c
	tmp_reaction = rate_constant_array(102)*(E_R_asn_deg)*(M_asn_L_c/(saturation_constant_array(102,122) + M_asn_L_c))*(M_amp_c/(saturation_constant_array(102,108) + M_amp_c))*(M_ppi_c/(saturation_constant_array(102,30) + M_ppi_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gly_deg: M_gly_L_c+M_accoa_c+M_o2_c+M_h2o_c = M_coa_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c
	tmp_reaction = rate_constant_array(103)*(E_R_gly_deg)*(M_gly_L_c/(saturation_constant_array(103,120) + M_gly_L_c))*(M_accoa_c/(saturation_constant_array(103,11) + M_accoa_c))*(M_o2_c/(saturation_constant_array(103,138) + M_o2_c))*(M_h2o_c/(saturation_constant_array(103,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mglx_deg: M_mglx_c+M_nad_c+M_h2o_c = M_pyr_c+M_nadh_c+M_h_c
	tmp_reaction = rate_constant_array(104)*(E_R_mglx_deg)*(M_mglx_c/(saturation_constant_array(104,44) + M_mglx_c))*(M_nad_c/(saturation_constant_array(104,145) + M_nad_c))*(M_h2o_c/(saturation_constant_array(104,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ser_deg: M_ser_L_c = M_nh3_c+M_pyr_c
	tmp_reaction = rate_constant_array(105)*(E_R_ser_deg)*(M_ser_L_c/(saturation_constant_array(105,129) + M_ser_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_pro_deg: M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c = M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
	tmp_reaction = rate_constant_array(106)*(E_R_pro_deg)*(M_pro_L_c/(saturation_constant_array(106,128) + M_pro_L_c))*(M_q8_c/(saturation_constant_array(106,24) + M_q8_c))*(M_h2o_c/(saturation_constant_array(106,144) + M_h2o_c))*(M_nad_c/(saturation_constant_array(106,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_thr_deg1: M_thr_L_c+M_nad_c+M_coa_c = M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
	tmp_reaction = rate_constant_array(107)*(E_R_thr_deg1)*(M_thr_L_c/(saturation_constant_array(107,130) + M_thr_L_c))*(M_nad_c/(saturation_constant_array(107,145) + M_nad_c))*(M_coa_c/(saturation_constant_array(107,10) + M_coa_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_thr_deg2: M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c = M_nadh_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c+M_h_c
	tmp_reaction = rate_constant_array(108)*(E_R_thr_deg2)*(M_thr_L_c/(saturation_constant_array(108,130) + M_thr_L_c))*(M_nad_c/(saturation_constant_array(108,145) + M_nad_c))*(M_o2_c/(saturation_constant_array(108,138) + M_o2_c))*(M_h2o_c/(saturation_constant_array(108,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_thr_deg3: M_thr_L_c+M_pi_c+M_adp_c = M_nh3_c+M_for_c+M_atp_c+M_prop_c
	tmp_reaction = rate_constant_array(109)*(E_R_thr_deg3)*(M_thr_L_c/(saturation_constant_array(109,130) + M_thr_L_c))*(M_pi_c/(saturation_constant_array(109,140) + M_pi_c))*(M_adp_c/(saturation_constant_array(109,107) + M_adp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_trp_deg: M_trp_L_c+M_h2o_c = M_indole_c+M_nh3_c+M_pyr_c
	tmp_reaction = rate_constant_array(110)*(E_R_trp_deg)*(M_trp_L_c/(saturation_constant_array(110,131) + M_trp_L_c))*(M_h2o_c/(saturation_constant_array(110,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_cys_deg: M_cys_L_c+M_h2o_c = M_h2s_c+M_nh3_c+M_pyr_c
	tmp_reaction = rate_constant_array(111)*(E_R_cys_deg)*(M_cys_L_c/(saturation_constant_array(111,123) + M_cys_L_c))*(M_h2o_c/(saturation_constant_array(111,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_lys_deg: M_lys_L_c = M_co2_c+M_cadav_c
	tmp_reaction = rate_constant_array(112)*(E_R_lys_deg)*(M_lys_L_c/(saturation_constant_array(112,124) + M_lys_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gln_deg: M_gln_L_c+M_h2o_c = M_nh3_c+M_glu_L_c
	tmp_reaction = rate_constant_array(113)*(E_R_gln_deg)*(M_gln_L_c/(saturation_constant_array(113,137) + M_gln_L_c))*(M_h2o_c/(saturation_constant_array(113,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_glu_deg: M_glu_L_c = M_co2_c+M_gaba_c
	tmp_reaction = rate_constant_array(114)*(E_R_glu_deg)*(M_glu_L_c/(saturation_constant_array(114,136) + M_glu_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gaba_deg1: M_gaba_c+M_akg_c+M_h2o_c+M_nad_c = M_succ_c+M_glu_L_c+M_h_c+M_nadh_c
	tmp_reaction = rate_constant_array(115)*(E_R_gaba_deg1)*(M_gaba_c/(saturation_constant_array(115,48) + M_gaba_c))*(M_akg_c/(saturation_constant_array(115,22) + M_akg_c))*(M_h2o_c/(saturation_constant_array(115,144) + M_h2o_c))*(M_nad_c/(saturation_constant_array(115,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gaba_deg2: M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c = M_succ_c+M_glu_L_c+M_h_c+M_nadph_c
	tmp_reaction = rate_constant_array(116)*(E_R_gaba_deg2)*(M_gaba_c/(saturation_constant_array(116,48) + M_gaba_c))*(M_akg_c/(saturation_constant_array(116,22) + M_akg_c))*(M_h2o_c/(saturation_constant_array(116,144) + M_h2o_c))*(M_nadp_c/(saturation_constant_array(116,147) + M_nadp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_chor: M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c+M_h_c = M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c
	tmp_reaction = rate_constant_array(117)*(E_R_chor)*(M_e4p_c/(saturation_constant_array(117,18) + M_e4p_c))*(M_pep_c/(saturation_constant_array(117,101) + M_pep_c))*(M_nadph_c/(saturation_constant_array(117,148) + M_nadph_c))*(M_atp_c/(saturation_constant_array(117,106) + M_atp_c))*(M_h_c/(saturation_constant_array(117,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fol_e: M_gtp_c+4*M_h2o_c = M_for_c+3.0*M_pi_c+M_glycoA_c+M_78mdp_c
	tmp_reaction = rate_constant_array(118)*(E_R_fol_e)*(M_gtp_c/(saturation_constant_array(118,109) + M_gtp_c))*(M_h2o_c/(saturation_constant_array(118,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fol_1: M_chor_c+M_gln_L_c = M_4adochor_c+M_glu_L_c
	tmp_reaction = rate_constant_array(119)*(E_R_fol_1)*(M_chor_c/(saturation_constant_array(119,42) + M_chor_c))*(M_gln_L_c/(saturation_constant_array(119,137) + M_gln_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fol_2a: M_4adochor_c = M_4abz_c+M_pyr_c
	tmp_reaction = rate_constant_array(120)*(E_R_fol_2a)*(M_4adochor_c/(saturation_constant_array(120,51) + M_4adochor_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fol_2b: M_4abz_c+M_78mdp_c = M_78dhf_c+M_h2o_c
	tmp_reaction = rate_constant_array(121)*(E_R_fol_2b)*(M_4abz_c/(saturation_constant_array(121,52) + M_4abz_c))*(M_78mdp_c/(saturation_constant_array(121,50) + M_78mdp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fol_3: M_78dhf_c+M_atp_c+M_glu_L_c = M_adp_c+M_pi_c+M_dhf_c
	tmp_reaction = rate_constant_array(122)*(E_R_fol_3)*(M_78dhf_c/(saturation_constant_array(122,53) + M_78dhf_c))*(M_atp_c/(saturation_constant_array(122,106) + M_atp_c))*(M_glu_L_c/(saturation_constant_array(122,136) + M_glu_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_fol_4: M_dhf_c+M_nadph_c+M_h_c = M_thf_c+M_nadp_c
	tmp_reaction = rate_constant_array(123)*(E_R_fol_4)*(M_dhf_c/(saturation_constant_array(123,54) + M_dhf_c))*(M_nadph_c/(saturation_constant_array(123,148) + M_nadph_c))*(M_h_c/(saturation_constant_array(123,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gly_fol: M_gly_L_c+M_thf_c+M_nad_c = M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c
	tmp_reaction = rate_constant_array(124)*(E_R_gly_fol)*(M_gly_L_c/(saturation_constant_array(124,120) + M_gly_L_c))*(M_thf_c/(saturation_constant_array(124,38) + M_thf_c))*(M_nad_c/(saturation_constant_array(124,145) + M_nad_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gly_fol_R: M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c = M_gly_L_c+M_thf_c+M_nad_c
	tmp_reaction = rate_constant_array(125)*(E_R_gly_fol_R)*(M_mlthf_c/(saturation_constant_array(125,39) + M_mlthf_c))*(M_nh3_c/(saturation_constant_array(125,35) + M_nh3_c))*(M_co2_c/(saturation_constant_array(125,139) + M_co2_c))*(M_nadh_c/(saturation_constant_array(125,146) + M_nadh_c))*(M_h_c/(saturation_constant_array(125,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mthfd: M_mlthf_c+M_nadp_c = M_methf_c+M_nadph_c
	tmp_reaction = rate_constant_array(126)*(E_R_mthfd)*(M_mlthf_c/(saturation_constant_array(126,39) + M_mlthf_c))*(M_nadp_c/(saturation_constant_array(126,147) + M_nadp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mthfd_R: M_methf_c+M_nadph_c = M_mlthf_c+M_nadp_c
	tmp_reaction = rate_constant_array(127)*(E_R_mthfd_R)*(M_methf_c/(saturation_constant_array(127,55) + M_methf_c))*(M_nadph_c/(saturation_constant_array(127,148) + M_nadph_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mthfc: M_h2o_c+M_methf_c = M_10fthf_c+M_h_c
	tmp_reaction = rate_constant_array(128)*(E_R_mthfc)*(M_h2o_c/(saturation_constant_array(128,144) + M_h2o_c))*(M_methf_c/(saturation_constant_array(128,55) + M_methf_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mthfc_R: M_10fthf_c+M_h_c = M_h2o_c+M_methf_c
	tmp_reaction = rate_constant_array(129)*(E_R_mthfc_R)*(M_10fthf_c/(saturation_constant_array(129,56) + M_10fthf_c))*(M_h_c/(saturation_constant_array(129,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mthfr2a: M_mlthf_c+M_h_c+M_nadh_c = M_5mthf_c+M_nad_c
	tmp_reaction = rate_constant_array(130)*(E_R_mthfr2a)*(M_mlthf_c/(saturation_constant_array(130,39) + M_mlthf_c))*(M_h_c/(saturation_constant_array(130,143) + M_h_c))*(M_nadh_c/(saturation_constant_array(130,146) + M_nadh_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_mthfr2b: M_mlthf_c+M_h_c+M_nadph_c = M_5mthf_c+M_nadp_c
	tmp_reaction = rate_constant_array(131)*(E_R_mthfr2b)*(M_mlthf_c/(saturation_constant_array(131,39) + M_mlthf_c))*(M_h_c/(saturation_constant_array(131,143) + M_h_c))*(M_nadph_c/(saturation_constant_array(131,148) + M_nadph_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_prpp_syn: M_r5p_c+M_atp_c = M_prpp_c+M_amp_c
	tmp_reaction = rate_constant_array(132)*(E_R_prpp_syn)*(M_r5p_c/(saturation_constant_array(132,16) + M_r5p_c))*(M_atp_c/(saturation_constant_array(132,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_or_syn_1: 2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c+M_h_c = 2.0*M_adp_c+M_glu_L_c+M_pi_c+M_clasp_c
	tmp_reaction = rate_constant_array(133)*(E_R_or_syn_1)*(M_atp_c/(saturation_constant_array(133,106) + M_atp_c))*(M_gln_L_c/(saturation_constant_array(133,137) + M_gln_L_c))*(M_hco3_c/(saturation_constant_array(133,58) + M_hco3_c))*(M_h2o_c/(saturation_constant_array(133,144) + M_h2o_c))*(M_h_c/(saturation_constant_array(133,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_or_syn_2: M_clasp_c+M_asp_L_c+M_q8_c = M_or_c+M_q8h2_c+M_h2o_c+M_pi_c
	tmp_reaction = rate_constant_array(134)*(E_R_or_syn_2)*(M_clasp_c/(saturation_constant_array(134,59) + M_clasp_c))*(M_asp_L_c/(saturation_constant_array(134,119) + M_asp_L_c))*(M_q8_c/(saturation_constant_array(134,24) + M_q8_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_omp_syn: M_prpp_c+M_or_c = M_omp_c+M_ppi_c
	tmp_reaction = rate_constant_array(135)*(E_R_omp_syn)*(M_prpp_c/(saturation_constant_array(135,57) + M_prpp_c))*(M_or_c/(saturation_constant_array(135,60) + M_or_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ump_syn: M_omp_c = M_ump_c+M_co2_c
	tmp_reaction = rate_constant_array(136)*(E_R_ump_syn)*(M_omp_c/(saturation_constant_array(136,61) + M_omp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ctp_1: M_utp_c+M_atp_c+M_nh3_c = M_ctp_c+M_adp_c+M_pi_c
	tmp_reaction = rate_constant_array(137)*(E_R_ctp_1)*(M_utp_c/(saturation_constant_array(137,112) + M_utp_c))*(M_atp_c/(saturation_constant_array(137,106) + M_atp_c))*(M_nh3_c/(saturation_constant_array(137,35) + M_nh3_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ctp_2: M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c = M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c
	tmp_reaction = rate_constant_array(138)*(E_R_ctp_2)*(M_utp_c/(saturation_constant_array(138,112) + M_utp_c))*(M_gln_L_c/(saturation_constant_array(138,137) + M_gln_L_c))*(M_atp_c/(saturation_constant_array(138,106) + M_atp_c))*(M_h2o_c/(saturation_constant_array(138,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_1: M_gln_L_c+M_prpp_c+M_h2o_c = M_5pbdra+M_ppi_c+M_glu_L_c
	tmp_reaction = rate_constant_array(139)*(E_R_A_syn_1)*(M_gln_L_c/(saturation_constant_array(139,137) + M_gln_L_c))*(M_prpp_c/(saturation_constant_array(139,57) + M_prpp_c))*(M_h2o_c/(saturation_constant_array(139,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_2: M_atp_c+M_5pbdra+M_gly_L_c = M_adp_c+M_pi_c+M_gar_c
	tmp_reaction = rate_constant_array(140)*(E_R_A_syn_2)*(M_atp_c/(saturation_constant_array(140,106) + M_atp_c))*(M_5pbdra/(saturation_constant_array(140,62) + M_5pbdra))*(M_gly_L_c/(saturation_constant_array(140,120) + M_gly_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_3: M_10fthf_c+M_gar_c = M_thf_c+M_fgar_c
	tmp_reaction = rate_constant_array(141)*(E_R_A_syn_3)*(M_10fthf_c/(saturation_constant_array(141,56) + M_10fthf_c))*(M_gar_c/(saturation_constant_array(141,63) + M_gar_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_4: M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c = M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c
	tmp_reaction = rate_constant_array(142)*(E_R_A_syn_4)*(M_atp_c/(saturation_constant_array(142,106) + M_atp_c))*(M_fgar_c/(saturation_constant_array(142,64) + M_fgar_c))*(M_gln_L_c/(saturation_constant_array(142,137) + M_gln_L_c))*(M_h2o_c/(saturation_constant_array(142,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_5: M_atp_c+M_fgam_c = M_adp_c+M_pi_c+M_air_c
	tmp_reaction = rate_constant_array(143)*(E_R_A_syn_5)*(M_atp_c/(saturation_constant_array(143,106) + M_atp_c))*(M_fgam_c/(saturation_constant_array(143,65) + M_fgam_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_6: M_atp_c+M_air_c+M_hco3_c+M_h_c = M_adp_c+M_pi_c+M_cair_c
	tmp_reaction = rate_constant_array(144)*(E_R_A_syn_6)*(M_atp_c/(saturation_constant_array(144,106) + M_atp_c))*(M_air_c/(saturation_constant_array(144,66) + M_air_c))*(M_hco3_c/(saturation_constant_array(144,58) + M_hco3_c))*(M_h_c/(saturation_constant_array(144,143) + M_h_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_7: M_atp_c+M_cair_c+M_asp_L_c = M_adp_c+M_pi_c+M_saicar_c
	tmp_reaction = rate_constant_array(145)*(E_R_A_syn_7)*(M_atp_c/(saturation_constant_array(145,106) + M_atp_c))*(M_cair_c/(saturation_constant_array(145,67) + M_cair_c))*(M_asp_L_c/(saturation_constant_array(145,119) + M_asp_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_8: M_saicar_c = M_fum_c+M_aicar_c
	tmp_reaction = rate_constant_array(146)*(E_R_A_syn_8)*(M_saicar_c/(saturation_constant_array(146,68) + M_saicar_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_9: M_aicar_c+M_10fthf_c = M_thf_c+M_faicar_c
	tmp_reaction = rate_constant_array(147)*(E_R_A_syn_9)*(M_aicar_c/(saturation_constant_array(147,40) + M_aicar_c))*(M_10fthf_c/(saturation_constant_array(147,56) + M_10fthf_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_10: M_faicar_c = M_imp_c+M_h2o_c
	tmp_reaction = rate_constant_array(148)*(E_R_A_syn_10)*(M_faicar_c/(saturation_constant_array(148,69) + M_faicar_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_A_syn_12: M_imp_c+M_gtp_c+M_asp_L_c = M_gdp_c+M_pi_c+M_fum_c+M_amp_c
	tmp_reaction = rate_constant_array(149)*(E_R_A_syn_12)*(M_imp_c/(saturation_constant_array(149,70) + M_imp_c))*(M_gtp_c/(saturation_constant_array(149,109) + M_gtp_c))*(M_asp_L_c/(saturation_constant_array(149,119) + M_asp_L_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_xmp_syn: M_imp_c+M_nad_c+M_h2o_c = M_xmp_c+M_nadh_c+M_h_c
	tmp_reaction = rate_constant_array(150)*(E_R_xmp_syn)*(M_imp_c/(saturation_constant_array(150,70) + M_imp_c))*(M_nad_c/(saturation_constant_array(150,145) + M_nad_c))*(M_h2o_c/(saturation_constant_array(150,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gmp_syn: M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c = M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c
	tmp_reaction = rate_constant_array(151)*(E_R_gmp_syn)*(M_atp_c/(saturation_constant_array(151,106) + M_atp_c))*(M_xmp_c/(saturation_constant_array(151,71) + M_xmp_c))*(M_gln_L_c/(saturation_constant_array(151,137) + M_gln_L_c))*(M_h2o_c/(saturation_constant_array(151,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_atp_amp: M_atp_c+M_h2o_c = M_amp_c+M_ppi_c
	tmp_reaction = rate_constant_array(152)*(E_R_atp_amp)*(M_atp_c/(saturation_constant_array(152,106) + M_atp_c))*(M_h2o_c/(saturation_constant_array(152,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_utp_ump: M_utp_c+M_h2o_c = M_ump_c+M_ppi_c
	tmp_reaction = rate_constant_array(153)*(E_R_utp_ump)*(M_utp_c/(saturation_constant_array(153,112) + M_utp_c))*(M_h2o_c/(saturation_constant_array(153,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ctp_cmp: M_ctp_c+M_h2o_c = M_cmp_c+M_ppi_c
	tmp_reaction = rate_constant_array(154)*(E_R_ctp_cmp)*(M_ctp_c/(saturation_constant_array(154,115) + M_ctp_c))*(M_h2o_c/(saturation_constant_array(154,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gtp_gmp: M_gtp_c+M_h2o_c = M_gmp_c+M_ppi_c
	tmp_reaction = rate_constant_array(155)*(E_R_gtp_gmp)*(M_gtp_c/(saturation_constant_array(155,109) + M_gtp_c))*(M_h2o_c/(saturation_constant_array(155,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_atp_adp: M_atp_c+M_h2o_c = M_adp_c+M_pi_c
	tmp_reaction = rate_constant_array(156)*(E_R_atp_adp)*(M_atp_c/(saturation_constant_array(156,106) + M_atp_c))*(M_h2o_c/(saturation_constant_array(156,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_utp_adp: M_utp_c+M_h2o_c = M_udp_c+M_pi_c
	tmp_reaction = rate_constant_array(157)*(E_R_utp_adp)*(M_utp_c/(saturation_constant_array(157,112) + M_utp_c))*(M_h2o_c/(saturation_constant_array(157,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_ctp_adp: M_ctp_c+M_h2o_c = M_cdp_c+M_pi_c
	tmp_reaction = rate_constant_array(158)*(E_R_ctp_adp)*(M_ctp_c/(saturation_constant_array(158,115) + M_ctp_c))*(M_h2o_c/(saturation_constant_array(158,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gtp_adp: M_gtp_c+M_h2o_c = M_gdp_c+M_pi_c
	tmp_reaction = rate_constant_array(159)*(E_R_gtp_adp)*(M_gtp_c/(saturation_constant_array(159,109) + M_gtp_c))*(M_h2o_c/(saturation_constant_array(159,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_udp_utp: M_udp_c+M_atp_c = M_utp_c+M_adp_c
	tmp_reaction = rate_constant_array(160)*(E_R_udp_utp)*(M_udp_c/(saturation_constant_array(160,113) + M_udp_c))*(M_atp_c/(saturation_constant_array(160,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_cdp_ctp: M_cdp_c+M_atp_c = M_ctp_c+M_adp_c
	tmp_reaction = rate_constant_array(161)*(E_R_cdp_ctp)*(M_cdp_c/(saturation_constant_array(161,116) + M_cdp_c))*(M_atp_c/(saturation_constant_array(161,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_gdp_gtp: M_gdp_c+M_atp_c = M_gtp_c+M_adp_c
	tmp_reaction = rate_constant_array(162)*(E_R_gdp_gtp)*(M_gdp_c/(saturation_constant_array(162,110) + M_gdp_c))*(M_atp_c/(saturation_constant_array(162,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_atp_ump: M_atp_c+M_ump_c = M_adp_c+M_udp_c
	tmp_reaction = rate_constant_array(163)*(E_R_atp_ump)*(M_atp_c/(saturation_constant_array(163,106) + M_atp_c))*(M_ump_c/(saturation_constant_array(163,114) + M_ump_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_atp_cmp: M_atp_c+M_cmp_c = M_adp_c+M_cdp_c
	tmp_reaction = rate_constant_array(164)*(E_R_atp_cmp)*(M_atp_c/(saturation_constant_array(164,106) + M_atp_c))*(M_cmp_c/(saturation_constant_array(164,117) + M_cmp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_atp_gmp: M_atp_c+M_gmp_c = M_adp_c+M_gdp_c
	tmp_reaction = rate_constant_array(165)*(E_R_atp_gmp)*(M_atp_c/(saturation_constant_array(165,106) + M_atp_c))*(M_gmp_c/(saturation_constant_array(165,111) + M_gmp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% R_adk_atp: M_amp_c+M_atp_c = 2.0*M_adp_c
	tmp_reaction = rate_constant_array(166)*(E_R_adk_atp)*(M_amp_c/(saturation_constant_array(166,108) + M_amp_c))*(M_atp_c/(saturation_constant_array(166,106) + M_atp_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Import_o2: () = M_o2_c
	tmp_reaction = rate_constant_array(167)*(E_Import_o2);
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Import_co2: () = M_co2_c
	tmp_reaction = rate_constant_array(168)*(E_Import_co2);
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Import_pi: () = M_pi_c
	tmp_reaction = rate_constant_array(169)*(E_Import_pi);
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Import_nh4: () = M_nh4_c
	tmp_reaction = rate_constant_array(170)*(E_Import_nh4);
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Import_so4: () = M_so4_c
	tmp_reaction = rate_constant_array(171)*(E_Import_so4);
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Import_h2o: () = M_h2o_c
	tmp_reaction = rate_constant_array(172)*(E_Import_h2o);
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Export_o2: M_o2_c = ()
	tmp_reaction = rate_constant_array(173)*(E_Export_o2)*(M_o2_c/(saturation_constant_array(173,138) + M_o2_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Export_co2: M_co2_c = ()
	tmp_reaction = rate_constant_array(174)*(E_Export_co2)*(M_co2_c/(saturation_constant_array(174,139) + M_co2_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Export_pi: M_pi_c = ()
	tmp_reaction = rate_constant_array(175)*(E_Export_pi)*(M_pi_c/(saturation_constant_array(175,140) + M_pi_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Export_nh4: M_nh4_c = ()
	tmp_reaction = rate_constant_array(176)*(E_Export_nh4)*(M_nh4_c/(saturation_constant_array(176,141) + M_nh4_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Export_so4: M_so4_c = ()
	tmp_reaction = rate_constant_array(177)*(E_Export_so4)*(M_so4_c/(saturation_constant_array(177,142) + M_so4_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Export_h2o: M_h2o_c = ()
	tmp_reaction = rate_constant_array(178)*(E_Export_h2o)*(M_h2o_c/(saturation_constant_array(178,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Proton_gradient: M_h_e = M_h_c
	tmp_reaction = rate_constant_array(179)*(E_Proton_gradient)*(M_h_e/(saturation_constant_array(179,29) + M_h_e));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% transcriptional_initiation_CAT: GENE_CAT+RNAP = OPEN_GENE_CAT
	tmp_reaction = 0;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% transcription_CAT: GENE_CAT+RNAP+151*M_gtp_c+144*M_ctp_c+189*M_utp_c+176*M_atp_c+660*M_h2o_c = mRNA_CAT+GENE_CAT+RNAP+660*M_ppi_c
	tmp_reaction = rate_constant_array(181)*(E_transcription_CAT)*(GENE_CAT/(saturation_constant_array(180,72) + GENE_CAT))*(RNAP)*(M_gtp_c/(saturation_constant_array(181,109) + M_gtp_c))*(M_ctp_c/(saturation_constant_array(181,115) + M_ctp_c))*(M_utp_c/(saturation_constant_array(181,112) + M_utp_c))*(M_atp_c/(saturation_constant_array(181,106) + M_atp_c))*(M_h2o_c/(saturation_constant_array(181,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% mRNA_degradation_CAT: mRNA_CAT = 151*M_gmp_c+144*M_cmp_c+189*M_ump_c+176*M_amp_c
	tmp_reaction = rate_constant_array(182)*(E_mRNA_degradation_CAT)*(mRNA_CAT);
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% translation_initiation_CAT: mRNA_CAT+RIBOSOME = RIBOSOME_START_CAT
	tmp_reaction = 0;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% translation_CAT: mRNA_CAT+RIBOSOME+438*M_gtp_c+438*M_h2o_c+15.0*M_ala_L_c_tRNA+5.0*M_arg_L_c_tRNA+10.0*M_asn_L_c_tRNA+12.0*M_asp_L_c_tRNA+5.0*M_cys_L_c_tRNA+12.0*M_glu_L_c_tRNA+13.0*M_gln_L_c_tRNA+10.0*M_gly_L_c_tRNA+12.0*M_his_L_c_tRNA+9.0*M_ile_L_c_tRNA+13.0*M_leu_L_c_tRNA+12.0*M_lys_L_c_tRNA+9.0*M_met_L_c_tRNA+20.0*M_phe_L_c_tRNA+7.0*M_pro_L_c_tRNA+10.0*M_ser_L_c_tRNA+13.0*M_thr_L_c_tRNA+5.0*M_trp_L_c_tRNA+11.0*M_tyr_L_c_tRNA+16.0*M_val_L_c_tRNA = RIBOSOME+mRNA_CAT+PROTEIN_CAT+438*M_gdp_c+438*M_pi_c+219*tRNA
	tmp_reaction = rate_constant_array(184)*(E_translation_CAT)*(mRNA_CAT/(saturation_constant_array(183,75) + mRNA_CAT))*(RIBOSOME)*(M_gtp_c/(saturation_constant_array(184,109) + M_gtp_c))*(M_h2o_c/(saturation_constant_array(184,144) + M_h2o_c))*(M_ala_L_c_tRNA/(saturation_constant_array(184,78) + M_ala_L_c_tRNA))*(M_arg_L_c_tRNA/(saturation_constant_array(184,79) + M_arg_L_c_tRNA))*(M_asn_L_c_tRNA/(saturation_constant_array(184,80) + M_asn_L_c_tRNA))*(M_asp_L_c_tRNA/(saturation_constant_array(184,81) + M_asp_L_c_tRNA))*(M_cys_L_c_tRNA/(saturation_constant_array(184,82) + M_cys_L_c_tRNA))*(M_glu_L_c_tRNA/(saturation_constant_array(184,83) + M_glu_L_c_tRNA))*(M_gln_L_c_tRNA/(saturation_constant_array(184,84) + M_gln_L_c_tRNA))*(M_gly_L_c_tRNA/(saturation_constant_array(184,85) + M_gly_L_c_tRNA))*(M_his_L_c_tRNA/(saturation_constant_array(184,86) + M_his_L_c_tRNA))*(M_ile_L_c_tRNA/(saturation_constant_array(184,87) + M_ile_L_c_tRNA))*(M_leu_L_c_tRNA/(saturation_constant_array(184,88) + M_leu_L_c_tRNA))*(M_lys_L_c_tRNA/(saturation_constant_array(184,89) + M_lys_L_c_tRNA))*(M_met_L_c_tRNA/(saturation_constant_array(184,90) + M_met_L_c_tRNA))*(M_phe_L_c_tRNA/(saturation_constant_array(184,91) + M_phe_L_c_tRNA))*(M_pro_L_c_tRNA/(saturation_constant_array(184,92) + M_pro_L_c_tRNA))*(M_ser_L_c_tRNA/(saturation_constant_array(184,93) + M_ser_L_c_tRNA))*(M_thr_L_c_tRNA/(saturation_constant_array(184,94) + M_thr_L_c_tRNA))*(M_trp_L_c_tRNA/(saturation_constant_array(184,95) + M_trp_L_c_tRNA))*(M_tyr_L_c_tRNA/(saturation_constant_array(184,96) + M_tyr_L_c_tRNA))*(M_val_L_c_tRNA/(saturation_constant_array(184,97) + M_val_L_c_tRNA));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_ala_L_c_CAT: 15.0*M_ala_L_c+15.0*M_atp_c+15.0*tRNA+15*M_h2o_c = 15.0*M_ala_L_c_tRNA+15.0*M_amp_c+15.0*M_ppi_c
	tmp_reaction = rate_constant_array(185)*(E_tRNA_charging_M_ala_L_c_CAT)*(M_ala_L_c/(saturation_constant_array(185,126) + M_ala_L_c))*(M_atp_c/(saturation_constant_array(185,106) + M_atp_c))*(tRNA/(saturation_constant_array(185,99) + tRNA))*(M_h2o_c/(saturation_constant_array(185,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_arg_L_c_CAT: 5.0*M_arg_L_c+5.0*M_atp_c+5.0*tRNA+5*M_h2o_c = 5.0*M_arg_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
	tmp_reaction = rate_constant_array(186)*(E_tRNA_charging_M_arg_L_c_CAT)*(M_arg_L_c/(saturation_constant_array(186,36) + M_arg_L_c))*(M_atp_c/(saturation_constant_array(186,106) + M_atp_c))*(tRNA/(saturation_constant_array(186,99) + tRNA))*(M_h2o_c/(saturation_constant_array(186,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_asn_L_c_CAT: 10.0*M_asn_L_c+10.0*M_atp_c+10.0*tRNA+10*M_h2o_c = 10.0*M_asn_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
	tmp_reaction = rate_constant_array(187)*(E_tRNA_charging_M_asn_L_c_CAT)*(M_asn_L_c/(saturation_constant_array(187,122) + M_asn_L_c))*(M_atp_c/(saturation_constant_array(187,106) + M_atp_c))*(tRNA/(saturation_constant_array(187,99) + tRNA))*(M_h2o_c/(saturation_constant_array(187,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_asp_L_c_CAT: 12.0*M_asp_L_c+12.0*M_atp_c+12.0*tRNA+12*M_h2o_c = 12.0*M_asp_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
	tmp_reaction = rate_constant_array(188)*(E_tRNA_charging_M_asp_L_c_CAT)*(M_asp_L_c/(saturation_constant_array(188,119) + M_asp_L_c))*(M_atp_c/(saturation_constant_array(188,106) + M_atp_c))*(tRNA/(saturation_constant_array(188,99) + tRNA))*(M_h2o_c/(saturation_constant_array(188,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_cys_L_c_CAT: 5.0*M_cys_L_c+5.0*M_atp_c+5.0*tRNA+5*M_h2o_c = 5.0*M_cys_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
	tmp_reaction = rate_constant_array(189)*(E_tRNA_charging_M_cys_L_c_CAT)*(M_cys_L_c/(saturation_constant_array(189,123) + M_cys_L_c))*(M_atp_c/(saturation_constant_array(189,106) + M_atp_c))*(tRNA/(saturation_constant_array(189,99) + tRNA))*(M_h2o_c/(saturation_constant_array(189,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_glu_L_c_CAT: 12.0*M_glu_L_c+12.0*M_atp_c+12.0*tRNA+12*M_h2o_c = 12.0*M_glu_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
	tmp_reaction = rate_constant_array(190)*(E_tRNA_charging_M_glu_L_c_CAT)*(M_glu_L_c/(saturation_constant_array(190,136) + M_glu_L_c))*(M_atp_c/(saturation_constant_array(190,106) + M_atp_c))*(tRNA/(saturation_constant_array(190,99) + tRNA))*(M_h2o_c/(saturation_constant_array(190,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_gln_L_c_CAT: 13.0*M_gln_L_c+13.0*M_atp_c+13.0*tRNA+13*M_h2o_c = 13.0*M_gln_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
	tmp_reaction = rate_constant_array(191)*(E_tRNA_charging_M_gln_L_c_CAT)*(M_gln_L_c/(saturation_constant_array(191,137) + M_gln_L_c))*(M_atp_c/(saturation_constant_array(191,106) + M_atp_c))*(tRNA/(saturation_constant_array(191,99) + tRNA))*(M_h2o_c/(saturation_constant_array(191,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_gly_L_c_CAT: 10.0*M_gly_L_c+10.0*M_atp_c+10.0*tRNA+10*M_h2o_c = 10.0*M_gly_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
	tmp_reaction = rate_constant_array(192)*(E_tRNA_charging_M_gly_L_c_CAT)*(M_gly_L_c/(saturation_constant_array(192,120) + M_gly_L_c))*(M_atp_c/(saturation_constant_array(192,106) + M_atp_c))*(tRNA/(saturation_constant_array(192,99) + tRNA))*(M_h2o_c/(saturation_constant_array(192,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_his_L_c_CAT: 12.0*M_his_L_c+12.0*M_atp_c+12.0*tRNA+12*M_h2o_c = 12.0*M_his_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
	tmp_reaction = rate_constant_array(193)*(E_tRNA_charging_M_his_L_c_CAT)*(M_his_L_c/(saturation_constant_array(193,125) + M_his_L_c))*(M_atp_c/(saturation_constant_array(193,106) + M_atp_c))*(tRNA/(saturation_constant_array(193,99) + tRNA))*(M_h2o_c/(saturation_constant_array(193,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_ile_L_c_CAT: 9.0*M_ile_L_c+9.0*M_atp_c+9.0*tRNA+9*M_h2o_c = 9.0*M_ile_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c
	tmp_reaction = rate_constant_array(194)*(E_tRNA_charging_M_ile_L_c_CAT)*(M_ile_L_c/(saturation_constant_array(194,121) + M_ile_L_c))*(M_atp_c/(saturation_constant_array(194,106) + M_atp_c))*(tRNA/(saturation_constant_array(194,99) + tRNA))*(M_h2o_c/(saturation_constant_array(194,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_leu_L_c_CAT: 13.0*M_leu_L_c+13.0*M_atp_c+13.0*tRNA+13*M_h2o_c = 13.0*M_leu_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
	tmp_reaction = rate_constant_array(195)*(E_tRNA_charging_M_leu_L_c_CAT)*(M_leu_L_c/(saturation_constant_array(195,135) + M_leu_L_c))*(M_atp_c/(saturation_constant_array(195,106) + M_atp_c))*(tRNA/(saturation_constant_array(195,99) + tRNA))*(M_h2o_c/(saturation_constant_array(195,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_lys_L_c_CAT: 12.0*M_lys_L_c+12.0*M_atp_c+12.0*tRNA+12*M_h2o_c = 12.0*M_lys_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
	tmp_reaction = rate_constant_array(196)*(E_tRNA_charging_M_lys_L_c_CAT)*(M_lys_L_c/(saturation_constant_array(196,124) + M_lys_L_c))*(M_atp_c/(saturation_constant_array(196,106) + M_atp_c))*(tRNA/(saturation_constant_array(196,99) + tRNA))*(M_h2o_c/(saturation_constant_array(196,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_met_L_c_CAT: 9.0*M_met_L_c+9.0*M_atp_c+9.0*tRNA+9*M_h2o_c = 9.0*M_met_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c
	tmp_reaction = rate_constant_array(197)*(E_tRNA_charging_M_met_L_c_CAT)*(M_met_L_c/(saturation_constant_array(197,134) + M_met_L_c))*(M_atp_c/(saturation_constant_array(197,106) + M_atp_c))*(tRNA/(saturation_constant_array(197,99) + tRNA))*(M_h2o_c/(saturation_constant_array(197,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_phe_L_c_CAT: 20.0*M_phe_L_c+20.0*M_atp_c+20.0*tRNA+20*M_h2o_c = 20.0*M_phe_L_c_tRNA+20.0*M_amp_c+20.0*M_ppi_c
	tmp_reaction = rate_constant_array(198)*(E_tRNA_charging_M_phe_L_c_CAT)*(M_phe_L_c/(saturation_constant_array(198,127) + M_phe_L_c))*(M_atp_c/(saturation_constant_array(198,106) + M_atp_c))*(tRNA/(saturation_constant_array(198,99) + tRNA))*(M_h2o_c/(saturation_constant_array(198,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_pro_L_c_CAT: 7.0*M_pro_L_c+7.0*M_atp_c+7.0*tRNA+7*M_h2o_c = 7.0*M_pro_L_c_tRNA+7.0*M_amp_c+7.0*M_ppi_c
	tmp_reaction = rate_constant_array(199)*(E_tRNA_charging_M_pro_L_c_CAT)*(M_pro_L_c/(saturation_constant_array(199,128) + M_pro_L_c))*(M_atp_c/(saturation_constant_array(199,106) + M_atp_c))*(tRNA/(saturation_constant_array(199,99) + tRNA))*(M_h2o_c/(saturation_constant_array(199,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_ser_L_c_CAT: 10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA+10*M_h2o_c = 10.0*M_ser_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
	tmp_reaction = rate_constant_array(200)*(E_tRNA_charging_M_ser_L_c_CAT)*(M_ser_L_c/(saturation_constant_array(200,129) + M_ser_L_c))*(M_atp_c/(saturation_constant_array(200,106) + M_atp_c))*(tRNA/(saturation_constant_array(200,99) + tRNA))*(M_h2o_c/(saturation_constant_array(200,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_thr_L_c_CAT: 13.0*M_thr_L_c+13.0*M_atp_c+13.0*tRNA+13*M_h2o_c = 13.0*M_thr_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
	tmp_reaction = rate_constant_array(201)*(E_tRNA_charging_M_thr_L_c_CAT)*(M_thr_L_c/(saturation_constant_array(201,130) + M_thr_L_c))*(M_atp_c/(saturation_constant_array(201,106) + M_atp_c))*(tRNA/(saturation_constant_array(201,99) + tRNA))*(M_h2o_c/(saturation_constant_array(201,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_trp_L_c_CAT: 5.0*M_trp_L_c+5.0*M_atp_c+5.0*tRNA+5*M_h2o_c = 5.0*M_trp_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
	tmp_reaction = rate_constant_array(202)*(E_tRNA_charging_M_trp_L_c_CAT)*(M_trp_L_c/(saturation_constant_array(202,131) + M_trp_L_c))*(M_atp_c/(saturation_constant_array(202,106) + M_atp_c))*(tRNA/(saturation_constant_array(202,99) + tRNA))*(M_h2o_c/(saturation_constant_array(202,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_tyr_L_c_CAT: 11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA+11*M_h2o_c = 11.0*M_tyr_L_c_tRNA+11.0*M_amp_c+11.0*M_ppi_c
	tmp_reaction = rate_constant_array(203)*(E_tRNA_charging_M_tyr_L_c_CAT)*(M_tyr_L_c/(saturation_constant_array(203,132) + M_tyr_L_c))*(M_atp_c/(saturation_constant_array(203,106) + M_atp_c))*(tRNA/(saturation_constant_array(203,99) + tRNA))*(M_h2o_c/(saturation_constant_array(203,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% tRNA_charging_M_val_L_c_CAT: 16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA+16*M_h2o_c = 16.0*M_val_L_c_tRNA+16.0*M_amp_c+16.0*M_ppi_c
	tmp_reaction = rate_constant_array(204)*(E_tRNA_charging_M_val_L_c_CAT)*(M_val_L_c/(saturation_constant_array(204,133) + M_val_L_c))*(M_atp_c/(saturation_constant_array(204,106) + M_atp_c))*(tRNA/(saturation_constant_array(204,99) + tRNA))*(M_h2o_c/(saturation_constant_array(204,144) + M_h2o_c));
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_glk_atp = ()
	tmp_reaction = rate_constant_array(205)*E_R_glk_atp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pgi = ()
	tmp_reaction = rate_constant_array(206)*E_R_pgi;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pgi_R = ()
	tmp_reaction = rate_constant_array(207)*E_R_pgi_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pfk = ()
	tmp_reaction = rate_constant_array(208)*E_R_pfk;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fdp = ()
	tmp_reaction = rate_constant_array(209)*E_R_fdp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fbaA = ()
	tmp_reaction = rate_constant_array(210)*E_R_fbaA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fbaA_R = ()
	tmp_reaction = rate_constant_array(211)*E_R_fbaA_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_tpiA = ()
	tmp_reaction = rate_constant_array(212)*E_R_tpiA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_tpiA_R = ()
	tmp_reaction = rate_constant_array(213)*E_R_tpiA_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gapA = ()
	tmp_reaction = rate_constant_array(214)*E_R_gapA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gapA_R = ()
	tmp_reaction = rate_constant_array(215)*E_R_gapA_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pgk = ()
	tmp_reaction = rate_constant_array(216)*E_R_pgk;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pgk_R = ()
	tmp_reaction = rate_constant_array(217)*E_R_pgk_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gpm = ()
	tmp_reaction = rate_constant_array(218)*E_R_gpm;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gpm_R = ()
	tmp_reaction = rate_constant_array(219)*E_R_gpm_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_eno = ()
	tmp_reaction = rate_constant_array(220)*E_R_eno;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_eno_R = ()
	tmp_reaction = rate_constant_array(221)*E_R_eno_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pyk = ()
	tmp_reaction = rate_constant_array(222)*E_R_pyk;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pck = ()
	tmp_reaction = rate_constant_array(223)*E_R_pck;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ppc = ()
	tmp_reaction = rate_constant_array(224)*E_R_ppc;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pdh = ()
	tmp_reaction = rate_constant_array(225)*E_R_pdh;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pps = ()
	tmp_reaction = rate_constant_array(226)*E_R_pps;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_zwf = ()
	tmp_reaction = rate_constant_array(227)*E_R_zwf;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_zwf_R = ()
	tmp_reaction = rate_constant_array(228)*E_R_zwf_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pgl = ()
	tmp_reaction = rate_constant_array(229)*E_R_pgl;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gnd = ()
	tmp_reaction = rate_constant_array(230)*E_R_gnd;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_rpe = ()
	tmp_reaction = rate_constant_array(231)*E_R_rpe;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_rpe_R = ()
	tmp_reaction = rate_constant_array(232)*E_R_rpe_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_rpi = ()
	tmp_reaction = rate_constant_array(233)*E_R_rpi;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_rpi_R = ()
	tmp_reaction = rate_constant_array(234)*E_R_rpi_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_talAB = ()
	tmp_reaction = rate_constant_array(235)*E_R_talAB;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_talAB_R = ()
	tmp_reaction = rate_constant_array(236)*E_R_talAB_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_tkt1 = ()
	tmp_reaction = rate_constant_array(237)*E_R_tkt1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_tkt1_R = ()
	tmp_reaction = rate_constant_array(238)*E_R_tkt1_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_tkt2 = ()
	tmp_reaction = rate_constant_array(239)*E_R_tkt2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_tkt2_R = ()
	tmp_reaction = rate_constant_array(240)*E_R_tkt2_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_edd = ()
	tmp_reaction = rate_constant_array(241)*E_R_edd;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_eda = ()
	tmp_reaction = rate_constant_array(242)*E_R_eda;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gltA = ()
	tmp_reaction = rate_constant_array(243)*E_R_gltA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_acn = ()
	tmp_reaction = rate_constant_array(244)*E_R_acn;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_acn_R = ()
	tmp_reaction = rate_constant_array(245)*E_R_acn_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_icd = ()
	tmp_reaction = rate_constant_array(246)*E_R_icd;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_icd_R = ()
	tmp_reaction = rate_constant_array(247)*E_R_icd_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_sucAB = ()
	tmp_reaction = rate_constant_array(248)*E_R_sucAB;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_sucCD = ()
	tmp_reaction = rate_constant_array(249)*E_R_sucCD;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_sdh = ()
	tmp_reaction = rate_constant_array(250)*E_R_sdh;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_frd = ()
	tmp_reaction = rate_constant_array(251)*E_R_frd;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fum = ()
	tmp_reaction = rate_constant_array(252)*E_R_fum;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fum_R = ()
	tmp_reaction = rate_constant_array(253)*E_R_fum_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mdh = ()
	tmp_reaction = rate_constant_array(254)*E_R_mdh;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mdh_R = ()
	tmp_reaction = rate_constant_array(255)*E_R_mdh_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_cyd = ()
	tmp_reaction = rate_constant_array(256)*E_R_cyd;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_cyo = ()
	tmp_reaction = rate_constant_array(257)*E_R_cyo;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_app = ()
	tmp_reaction = rate_constant_array(258)*E_R_app;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_atp = ()
	tmp_reaction = rate_constant_array(259)*E_R_atp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_nuo = ()
	tmp_reaction = rate_constant_array(260)*E_R_nuo;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pnt1 = ()
	tmp_reaction = rate_constant_array(261)*E_R_pnt1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pnt2 = ()
	tmp_reaction = rate_constant_array(262)*E_R_pnt2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ndh1 = ()
	tmp_reaction = rate_constant_array(263)*E_R_ndh1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ndh2 = ()
	tmp_reaction = rate_constant_array(264)*E_R_ndh2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ppa = ()
	tmp_reaction = rate_constant_array(265)*E_R_ppa;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_aceA = ()
	tmp_reaction = rate_constant_array(266)*E_R_aceA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_aceB = ()
	tmp_reaction = rate_constant_array(267)*E_R_aceB;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_maeA = ()
	tmp_reaction = rate_constant_array(268)*E_R_maeA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_maeB = ()
	tmp_reaction = rate_constant_array(269)*E_R_maeB;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pta = ()
	tmp_reaction = rate_constant_array(270)*E_R_pta;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pta_R = ()
	tmp_reaction = rate_constant_array(271)*E_R_pta_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ackA = ()
	tmp_reaction = rate_constant_array(272)*E_R_ackA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ackA_R = ()
	tmp_reaction = rate_constant_array(273)*E_R_ackA_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_acs = ()
	tmp_reaction = rate_constant_array(274)*E_R_acs;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_adhE = ()
	tmp_reaction = rate_constant_array(275)*E_R_adhE;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_adhE_R = ()
	tmp_reaction = rate_constant_array(276)*E_R_adhE_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ldh = ()
	tmp_reaction = rate_constant_array(277)*E_R_ldh;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ldh_R = ()
	tmp_reaction = rate_constant_array(278)*E_R_ldh_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pflAB = ()
	tmp_reaction = rate_constant_array(279)*E_R_pflAB;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_alaAC = ()
	tmp_reaction = rate_constant_array(280)*E_R_alaAC;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_alaAC_R = ()
	tmp_reaction = rate_constant_array(281)*E_R_alaAC_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_arg = ()
	tmp_reaction = rate_constant_array(282)*E_R_arg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_aspC = ()
	tmp_reaction = rate_constant_array(283)*E_R_aspC;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_asnB = ()
	tmp_reaction = rate_constant_array(284)*E_R_asnB;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_asnA = ()
	tmp_reaction = rate_constant_array(285)*E_R_asnA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_cysEMK = ()
	tmp_reaction = rate_constant_array(286)*E_R_cysEMK;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gltBD = ()
	tmp_reaction = rate_constant_array(287)*E_R_gltBD;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gdhA = ()
	tmp_reaction = rate_constant_array(288)*E_R_gdhA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gdhA_R = ()
	tmp_reaction = rate_constant_array(289)*E_R_gdhA_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_glnA = ()
	tmp_reaction = rate_constant_array(290)*E_R_glnA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_glyA = ()
	tmp_reaction = rate_constant_array(291)*E_R_glyA;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_his = ()
	tmp_reaction = rate_constant_array(292)*E_R_his;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ile = ()
	tmp_reaction = rate_constant_array(293)*E_R_ile;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_leu = ()
	tmp_reaction = rate_constant_array(294)*E_R_leu;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_lys = ()
	tmp_reaction = rate_constant_array(295)*E_R_lys;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_met = ()
	tmp_reaction = rate_constant_array(296)*E_R_met;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_phe = ()
	tmp_reaction = rate_constant_array(297)*E_R_phe;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pro = ()
	tmp_reaction = rate_constant_array(298)*E_R_pro;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_serABC = ()
	tmp_reaction = rate_constant_array(299)*E_R_serABC;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_thr = ()
	tmp_reaction = rate_constant_array(300)*E_R_thr;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_trp = ()
	tmp_reaction = rate_constant_array(301)*E_R_trp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_tyr = ()
	tmp_reaction = rate_constant_array(302)*E_R_tyr;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_val = ()
	tmp_reaction = rate_constant_array(303)*E_R_val;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_arg_deg = ()
	tmp_reaction = rate_constant_array(304)*E_R_arg_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_asp_deg = ()
	tmp_reaction = rate_constant_array(305)*E_R_asp_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_asn_deg = ()
	tmp_reaction = rate_constant_array(306)*E_R_asn_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gly_deg = ()
	tmp_reaction = rate_constant_array(307)*E_R_gly_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mglx_deg = ()
	tmp_reaction = rate_constant_array(308)*E_R_mglx_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ser_deg = ()
	tmp_reaction = rate_constant_array(309)*E_R_ser_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_pro_deg = ()
	tmp_reaction = rate_constant_array(310)*E_R_pro_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_thr_deg1 = ()
	tmp_reaction = rate_constant_array(311)*E_R_thr_deg1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_thr_deg2 = ()
	tmp_reaction = rate_constant_array(312)*E_R_thr_deg2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_thr_deg3 = ()
	tmp_reaction = rate_constant_array(313)*E_R_thr_deg3;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_trp_deg = ()
	tmp_reaction = rate_constant_array(314)*E_R_trp_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_cys_deg = ()
	tmp_reaction = rate_constant_array(315)*E_R_cys_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_lys_deg = ()
	tmp_reaction = rate_constant_array(316)*E_R_lys_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gln_deg = ()
	tmp_reaction = rate_constant_array(317)*E_R_gln_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_glu_deg = ()
	tmp_reaction = rate_constant_array(318)*E_R_glu_deg;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gaba_deg1 = ()
	tmp_reaction = rate_constant_array(319)*E_R_gaba_deg1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gaba_deg2 = ()
	tmp_reaction = rate_constant_array(320)*E_R_gaba_deg2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_chor = ()
	tmp_reaction = rate_constant_array(321)*E_R_chor;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fol_e = ()
	tmp_reaction = rate_constant_array(322)*E_R_fol_e;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fol_1 = ()
	tmp_reaction = rate_constant_array(323)*E_R_fol_1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fol_2a = ()
	tmp_reaction = rate_constant_array(324)*E_R_fol_2a;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fol_2b = ()
	tmp_reaction = rate_constant_array(325)*E_R_fol_2b;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fol_3 = ()
	tmp_reaction = rate_constant_array(326)*E_R_fol_3;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_fol_4 = ()
	tmp_reaction = rate_constant_array(327)*E_R_fol_4;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gly_fol = ()
	tmp_reaction = rate_constant_array(328)*E_R_gly_fol;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gly_fol_R = ()
	tmp_reaction = rate_constant_array(329)*E_R_gly_fol_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mthfd = ()
	tmp_reaction = rate_constant_array(330)*E_R_mthfd;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mthfd_R = ()
	tmp_reaction = rate_constant_array(331)*E_R_mthfd_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mthfc = ()
	tmp_reaction = rate_constant_array(332)*E_R_mthfc;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mthfc_R = ()
	tmp_reaction = rate_constant_array(333)*E_R_mthfc_R;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mthfr2a = ()
	tmp_reaction = rate_constant_array(334)*E_R_mthfr2a;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_mthfr2b = ()
	tmp_reaction = rate_constant_array(335)*E_R_mthfr2b;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_prpp_syn = ()
	tmp_reaction = rate_constant_array(336)*E_R_prpp_syn;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_or_syn_1 = ()
	tmp_reaction = rate_constant_array(337)*E_R_or_syn_1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_or_syn_2 = ()
	tmp_reaction = rate_constant_array(338)*E_R_or_syn_2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_omp_syn = ()
	tmp_reaction = rate_constant_array(339)*E_R_omp_syn;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ump_syn = ()
	tmp_reaction = rate_constant_array(340)*E_R_ump_syn;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ctp_1 = ()
	tmp_reaction = rate_constant_array(341)*E_R_ctp_1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ctp_2 = ()
	tmp_reaction = rate_constant_array(342)*E_R_ctp_2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_1 = ()
	tmp_reaction = rate_constant_array(343)*E_R_A_syn_1;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_2 = ()
	tmp_reaction = rate_constant_array(344)*E_R_A_syn_2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_3 = ()
	tmp_reaction = rate_constant_array(345)*E_R_A_syn_3;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_4 = ()
	tmp_reaction = rate_constant_array(346)*E_R_A_syn_4;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_5 = ()
	tmp_reaction = rate_constant_array(347)*E_R_A_syn_5;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_6 = ()
	tmp_reaction = rate_constant_array(348)*E_R_A_syn_6;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_7 = ()
	tmp_reaction = rate_constant_array(349)*E_R_A_syn_7;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_8 = ()
	tmp_reaction = rate_constant_array(350)*E_R_A_syn_8;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_9 = ()
	tmp_reaction = rate_constant_array(351)*E_R_A_syn_9;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_10 = ()
	tmp_reaction = rate_constant_array(352)*E_R_A_syn_10;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_A_syn_12 = ()
	tmp_reaction = rate_constant_array(353)*E_R_A_syn_12;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_xmp_syn = ()
	tmp_reaction = rate_constant_array(354)*E_R_xmp_syn;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gmp_syn = ()
	tmp_reaction = rate_constant_array(355)*E_R_gmp_syn;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_atp_amp = ()
	tmp_reaction = rate_constant_array(356)*E_R_atp_amp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_utp_ump = ()
	tmp_reaction = rate_constant_array(357)*E_R_utp_ump;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ctp_cmp = ()
	tmp_reaction = rate_constant_array(358)*E_R_ctp_cmp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gtp_gmp = ()
	tmp_reaction = rate_constant_array(359)*E_R_gtp_gmp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_atp_adp = ()
	tmp_reaction = rate_constant_array(360)*E_R_atp_adp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_utp_adp = ()
	tmp_reaction = rate_constant_array(361)*E_R_utp_adp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_ctp_adp = ()
	tmp_reaction = rate_constant_array(362)*E_R_ctp_adp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gtp_adp = ()
	tmp_reaction = rate_constant_array(363)*E_R_gtp_adp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_udp_utp = ()
	tmp_reaction = rate_constant_array(364)*E_R_udp_utp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_cdp_ctp = ()
	tmp_reaction = rate_constant_array(365)*E_R_cdp_ctp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_gdp_gtp = ()
	tmp_reaction = rate_constant_array(366)*E_R_gdp_gtp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_atp_ump = ()
	tmp_reaction = rate_constant_array(367)*E_R_atp_ump;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_atp_cmp = ()
	tmp_reaction = rate_constant_array(368)*E_R_atp_cmp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_atp_gmp = ()
	tmp_reaction = rate_constant_array(369)*E_R_atp_gmp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_R_adk_atp = ()
	tmp_reaction = rate_constant_array(370)*E_R_adk_atp;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Import_o2 = ()
	tmp_reaction = rate_constant_array(371)*E_Import_o2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Import_co2 = ()
	tmp_reaction = rate_constant_array(372)*E_Import_co2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Import_pi = ()
	tmp_reaction = rate_constant_array(373)*E_Import_pi;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Import_nh4 = ()
	tmp_reaction = rate_constant_array(374)*E_Import_nh4;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Import_so4 = ()
	tmp_reaction = rate_constant_array(375)*E_Import_so4;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Import_h2o = ()
	tmp_reaction = rate_constant_array(376)*E_Import_h2o;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Export_o2 = ()
	tmp_reaction = rate_constant_array(377)*E_Export_o2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Export_co2 = ()
	tmp_reaction = rate_constant_array(378)*E_Export_co2;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Export_pi = ()
	tmp_reaction = rate_constant_array(379)*E_Export_pi;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Export_nh4 = ()
	tmp_reaction = rate_constant_array(380)*E_Export_nh4;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Export_so4 = ()
	tmp_reaction = rate_constant_array(381)*E_Export_so4;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Export_h2o = ()
	tmp_reaction = rate_constant_array(382)*E_Export_h2o;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_Proton_gradient = ()
	tmp_reaction = rate_constant_array(383)*E_Proton_gradient;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_transcriptional_initiation_CAT = ()
	tmp_reaction = rate_constant_array(384)*E_transcriptional_initiation_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_transcription_CAT = ()
	tmp_reaction = rate_constant_array(385)*E_transcription_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_mRNA_degradation_CAT = ()
	tmp_reaction = rate_constant_array(386)*E_mRNA_degradation_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_translation_initiation_CAT = ()
	tmp_reaction = rate_constant_array(387)*E_translation_initiation_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_translation_CAT = ()
	tmp_reaction = rate_constant_array(388)*E_translation_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_ala_L_c_CAT = ()
	tmp_reaction = rate_constant_array(389)*E_tRNA_charging_M_ala_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_arg_L_c_CAT = ()
	tmp_reaction = rate_constant_array(390)*E_tRNA_charging_M_arg_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_asn_L_c_CAT = ()
	tmp_reaction = rate_constant_array(391)*E_tRNA_charging_M_asn_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_asp_L_c_CAT = ()
	tmp_reaction = rate_constant_array(392)*E_tRNA_charging_M_asp_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_cys_L_c_CAT = ()
	tmp_reaction = rate_constant_array(393)*E_tRNA_charging_M_cys_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_glu_L_c_CAT = ()
	tmp_reaction = rate_constant_array(394)*E_tRNA_charging_M_glu_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_gln_L_c_CAT = ()
	tmp_reaction = rate_constant_array(395)*E_tRNA_charging_M_gln_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_gly_L_c_CAT = ()
	tmp_reaction = rate_constant_array(396)*E_tRNA_charging_M_gly_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_his_L_c_CAT = ()
	tmp_reaction = rate_constant_array(397)*E_tRNA_charging_M_his_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_ile_L_c_CAT = ()
	tmp_reaction = rate_constant_array(398)*E_tRNA_charging_M_ile_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_leu_L_c_CAT = ()
	tmp_reaction = rate_constant_array(399)*E_tRNA_charging_M_leu_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_lys_L_c_CAT = ()
	tmp_reaction = rate_constant_array(400)*E_tRNA_charging_M_lys_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_met_L_c_CAT = ()
	tmp_reaction = rate_constant_array(401)*E_tRNA_charging_M_met_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_phe_L_c_CAT = ()
	tmp_reaction = rate_constant_array(402)*E_tRNA_charging_M_phe_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_pro_L_c_CAT = ()
	tmp_reaction = rate_constant_array(403)*E_tRNA_charging_M_pro_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_ser_L_c_CAT = ()
	tmp_reaction = rate_constant_array(404)*E_tRNA_charging_M_ser_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_thr_L_c_CAT = ()
	tmp_reaction = rate_constant_array(405)*E_tRNA_charging_M_thr_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_trp_L_c_CAT = ()
	tmp_reaction = rate_constant_array(406)*E_tRNA_charging_M_trp_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_tyr_L_c_CAT = ()
	tmp_reaction = rate_constant_array(407)*E_tRNA_charging_M_tyr_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
	

	% Degradation: E_tRNA_charging_M_val_L_c_CAT = ()
	tmp_reaction = rate_constant_array(408)*E_tRNA_charging_M_val_L_c_CAT;
	rate_vector = [rate_vector; tmp_reaction];
	tmp_reaction = 0;
    
    dx = model.S * rate_vector;

end