# -*- coding:utf-8 -*-
import re
import json
import docxtpl
import os
import json
import copy
import datetime
from itertools import chain
from pypinyin import pinyin, Style
import itertools
from collections import defaultdict
from itertools import groupby
import sys
import locale
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import argparse
import jinja2
from jinja2 import filters

# 筛选诊断相关变异
def filter_diag(var_list):
    """
    筛选诊断相关变异
    :param var_list:
    :return:
    """
    if not var_list:
        return []
    diag_var = list(filter(lambda x: x["Diagnostic_str"] != '-', var_list))
    return diag_var
jinja2.filters.FILTERS["filter_diag"] = filter_diag

# 筛选治疗相关变异
def filter_pred(var_list):
    if not var_list:
        return []
    pred_var = list(filter(lambda x: x["regimen"] != [], var_list))
    return pred_var
jinja2.filters.FILTERS["filter_pred"] = filter_pred

# 转化字典-变异类型
type_strans = {"CorruptedStart": "起始密码子变异",
               "Extension": "延伸突变",
               "FrameShift_Deletion": "移码突变",
               "FrameShift_Duplication": "移码突变",
               "FrameShift_Insertion": "移码突变",
               "FrameShift_Substitution": "移码突变",
               "nonFrameShift_Deletion": "框内非移码突变",
               "nonFrameShift_Duplication": "框内非移码突变",
               "nonFrameShift_Insertion": "框内非移码突变",
               "nonFrameShift_Substitution": "非移码突变",
               "Nonsense_Mutation": "无义突变",
               "nonSynonymous_Substitution": "错义突变",
               "Splicing": "经典剪接位点变异",
               "Synonymous_Substitution": "同义突变",
               "5'UTR": "5'UTR区突变",
               "3'UTR": "3'UTR区突变",
               "FlankingRegion3": "侧翼区",
               "FlankingRegion5": "侧翼区",
               "Intronic": "内含子区突变"}
# 转化字典-胚系变异临床意义
clinic_strans_g = {
    "Pathogenic": 5,
    "Likely pathogenic": 4,
    "Uncertain": 3,
    "Likely benign": 2,
    "Benign": 1,
    "Oncogenic": 5,
    "Likely oncogenic": 4
}
# 转化字典-体细胞变异临床意义
clinic_strans_s = {
    "Pathogenic": 5,
    "Likely pathogenic": 4,
    "Uncertain": 3,
    "Likely benign": 2,
    "Benign": 1,
    "Oncogenic": 5,
    "Likely oncogenic": 4
}
# 转化字典-药物敏感性
sensi_stran = {
    "Sensitive": "敏感",
    "Resistant": "耐药",
    "Adverse Response": "不良反应",
    "Predicted-Sensitive": "预测敏感",
    "Predicted-Resistant": "预测耐药"
}
# 转化字典-将体细胞和胚系变异的临床意义转化为数字，便于变异排序
clinic_level = {
    "Pathogenic": 5,
    "Likely pathogenic": 4,
    "Uncertain": 3,
    "Likely benign": 2,
    "Benign": 1,
    "Oncogenic": 5,
    "Likely oncogenic": 4
}
# CNV基因对应转录本号（目前生信流程未返回该内容，先加在脚本里，后续更新）-v1.3.7
'''
cnv_transcript_cp40 = {
	"CDK4" : "NM_000075.4",
	"ERBB2" : "NM_004448.3",
	"MET" : "NM_000245.4",
	"MYC" : "NM_002467.6",
	"NKX2-1" : "NM_001079668.3",
	"EGFR" : "NM_005228.5"
}
'''
# 转录本更新 2024年4月11日 嵇梦晨
cnv_transcript_cp40 = {
    "CDK4" : "NM_000075.4",
	"ERBB2" : "NM_004448.4",
	"MET" : "NM_000245.4",
	"MYC" : "NM_002467.6",
	"NKX2-1" : "NM_001079668.3",
	"EGFR" : "NM_005228.5",
    #2026年1月15日，模板升级新增检测CNV基因转录本，MDM2, GLI1, AKT2, CCNE1, CDH1, CDK6, ESR1, FGF19, FGFR1, FGFR2, FGFR3, FLT3, GLI2, PDCD1LG2, RAF1, RICTOR, VEGFA
    "MDM2" : "NM_002392.5",
    "GLI1" : "NM_005269.4"，
    "AKT2" : "NM_001626.4",
    "CCNE1" : "NM_001238.5",
    "CDH1" : "NM_004360.5",
    "CDK6" : "NM_001259.4",
    "ESR1" : "NM_000125.4",
    "FGF19" : "NM_005257.4",
    "FGFR1" : "NM_023110.4",
    "FGFR2" : "NM_000141.4",
    "FGFR3" : "NM_000142.4",
    "FLT3" : "NM_004119.4",
    "GLI2" : "NM_005270.4",
    "PDCD1LG2" : "NM_014140.4",
    "RAF1" : "NM_002880.4",
    "RICTOR" : "NM_152756.4",
    "VEGFA" : "NM_001025366.3"
}

# 只看致癌性的基因
CARCINOGENICITY_GENE_LIST = [
    "ARAF", "CTNNB1", "CYP19A1", "ERBB2", "FGFR1", "FGFR3", "FGFR2",
    "FGFR4", "IDH1", "IDH2", "KDR", "MYC", "NFE2L2", "NOTCH1", "NTRK2",
    "NTRK3", "DDR2", "MAPK1", "ROS1", "NKX2-1", "KEAP1", "SF3B1"
]
# 只看致病性的基因
FUNCTION_GENE_LIST = [
    "APC", "ATM", "ATR", "BARD1", "BLM", "BMPR1A", "BRCA1", "BRCA2", "RUNX1",
    "CDH1", "CDKN1B", "CDKN2A", "CEBPA", "CHEK1", "CTNNA1", "ENG", "FANCA",
    "FANCC", "FH", "MSH6", "HDAC2", "EPCAM", "SMAD4", "MAX", "MEN1", "MITF",
    "MLH1", "MRE11", "MSH2", "MSH3", "MUTYH", "NBN", "NF1", "NF2", "PMS2",
    "POLD1", "POLE", "PPP2R2A", "PRKAR1A", "PTCH1", "PTEN", "RAD51", "RAD51C",
    "RAD51B", "RAD51D", "RB1", "RECQL", "SDHA", "SDHB", "SDHC", "SDHD", "SMARCA4",
    "SMARCB1", "STK11", "TP53", "TSC1", "TSC2", "VHL", "WRN", "WT1", "AXIN2",
    "BAP1", "RAD54L", "RECQL4", "RAD50", "CHEK2", "PALLD", "DICER1", "MLH3",
    "SUFU", "CDK12", "RNF43", "FANCL", "FANCI", "TMEM127", "ELAC2", "CDC73",
    "PALB2", "BRIP1", "ABRAXAS1", "FLCN", "GEN1", "SDHAF2",
]

evi_level_dict = {
	"Clinical-phase I" : "临床试验",
	"Clinical-phase II" : "临床试验",
	"Clinical-phase III" : "临床试验",
	"Clinical-phase IV" : "临床试验",
	"Clinical-retrospective" : "临床试验",
	"Clinical-unknown phase" : "临床试验",
	"Case report" : "案例报道",
	"Preclinical-in vitro" : "临床前证据",
	"Preclinical-in vivo" : "临床前证据"
}

class Base:
    def __init__(self, json_name, output_div):
        """
        传入的json_name 不带后缀
        :param json_name:
        :param output_div:
        """
        if not os.path.exists(output_div):
            os.makedirs(output_div)
        self.json_name = json_name
        self.data_js = json.load(open(os.path.join(output_div, json_name + '.json'), 'r', encoding='utf-8'))
        self.output_div = output_div

        # 填充结果保存
        self.data_js['result'] = dict()
        # 个人信息
        self.data_js['result']['sample_info'] = dict()
        # 结果小结
        self.data_js['result']['var_count'] = dict()
        # 变异结果汇总
        self.data_js['result']['var'] = dict()
        # QC
        self.data_js['result']['qc'] = dict()
        # 参考文献
        self.data_js['result']['reference'] = list()

    def _sample_id(self):
        """
        样本id
        :return:
        """
        # 不含FP
        sample_id = self.data_js['sample_info'].get('sample_id', '-')
        if sample_id != '-':
            if sample_id[0] == 'F':  # 有F开头的 F221310F001D01
                sample_id = 'F' + sample_id[1:].split('FD')[0].split('P')[0].split('F')[0]
            else:
                sample_id = sample_id.split('FD')[0].split('P')[0].split('F')[0].split('B')[0].split('E')[0]

        return sample_id
    
    def origin(self, var):
        return 'G' if var.get('var_origin') == 'germline' else 'S'

    def _drug_organization(self):
        """
        药物批准机构
        :return: self.data_js['result']['drug_dict']
        """
        # 药物对应批准机构
        drug_dict = dict()
        if self.data_js["drug"]:
            for drug_item in self.data_js["drug"]:
                if drug_item["general_name_cn"]:
                    drug_dict[drug_item["general_name_cn"]] = "/".join(drug_item["approval_organization"]) if drug_item[
                        "approval_organization"] else "-"
                if drug_item["general_name_en"]:
                    drug_dict[drug_item["general_name_en"]] = "/".join(drug_item["approval_organization"]) if drug_item[
                        "approval_organization"] else "-"
        self.data_js['result']['drug_dict'] = drug_dict
        return drug_dict

    def _topinyin(self, instr):
        """
        中文转拼音，用于中文排序
        :return:
        """
        trans = "".join(chain.from_iterable(pinyin(instr, Style.TONE3)))
        return trans
    
    def merge_evi(self, datainfo):
        merge_result = []
        if datainfo:
            tmp_dict = {}
            for evi in datainfo:
                tmp_dict.setdefault(evi["evi_interpretation"], [])
                tmp_dict[evi["evi_interpretation"]].append({"regimen_name" : evi["regimen_name"]})
            
            for k, v  in tmp_dict.items():
                merge_result.append({"regimen_name" : "、".join([i["regimen_name"] for i in v]),
                "evi_interpretation" : k})
        return merge_result
    
    def merge_diag(self, datainfo):
        merge_result = []
        if datainfo:
            tmp_dict = {}
            for diag in datainfo:
                tmp_dict.setdefault(diag['tumor_name_cn'], [])
                tmp_dict[diag['tumor_name_cn']].append({'evi_interpretation' : diag['evi_interpretation']})
            for k, v in tmp_dict.items():
                merge_result.append({'evi_interpretation' : '\n'.join([i['evi_interpretation'] for i in v]), 'tumor_name_cn' : k})
        return merge_result
    
    def sv_info(self, var):
        return var["five_prime_gene"]+":"+var["five_prime_region"]+"-"+var["three_prime_gene"]+":"+var["three_prime_region"]
    
    def sv_combination(self, rnasv):
        if not rnasv:
            return []
        # 处理rnasv
        rnasv = sorted(rnasv, key= lambda x: (x["five_prime_gene"], x["three_prime_gene"]))
        for var in rnasv:
            var['sv_info'] = self.sv_info(var)
        sv_sum = groupby(rnasv, lambda x: x['sv_info'])
        sv_combination = []
        for i, j in sv_sum:
            sv_group = list(j)
            major_sv = sv_group[0]
            major_sv['rnasv_reads'] = sum([int(float(var['rnasv_reads'])) for var in sv_group])
            sv_combination.append(major_sv)
        return sv_combination
    
    def sense_rule(self, clin_sign):
        return "0" if re.search("Sensitive", clin_sign) else "1" if re.search("Resistant", clin_sign) else clin_sign

    def _var_list(self):
        """
        变异列表
        :return: self.data_js['result']['var_list']
        """
        if self.data_js['result'].get('var_list', ''):  # 防止重复运算
            return self.data_js['result']['var_list']
        
        var_list = list()
        jsonDict = self.data_js
        regimen_dict, regimen_adaptation = self.get_regimen_appr(self.data_js["therapeutic_regimen"])
        regimen_rule = ["埃克替尼","厄洛替尼","吉非替尼","阿法替尼","达可替尼","奥希替尼","阿美替尼","伏美替尼","厄洛替尼+雷莫西尤单抗","厄洛替尼+贝伐珠单抗"]

        jsonDict['rna_sv'] = self.sv_combination(jsonDict['rna_sv'])
        # var为列表，元素为字典
        # 找出I、II、III类变异
        var = list(filter(lambda k: {k["clinical_significance"], k["function_classification"]} &
                                    {"Uncertain", "Likely pathogenic", "Pathogenic", "Oncogenic", "Likely oncogenic"},
                          jsonDict["snvindel"] + jsonDict["rna_sv"] + jsonDict["cnv"]))
        
        if not var:
            self.data_js['result']['var_list'] = list()
            return self.data_js['result']['var_list']
        
        for i in var:
            # 丰度
            i['freq_num'] = i['freq'] if 'freq' in i.keys() else i['cn_mean'] if 'cn_mean' in i.keys() else i['rnasv_reads'] if 'rnasv_reads' in i.keys() else 0
            
            if "freq" in i.keys() and i["freq"]:
                if re.search("%", str(i["freq"])):
                    i["freq"] = i["freq"]
                else:
                    i["freq"] = "{:.2%}".format(i["freq"])
            elif "cn_mean" in i.keys() and i["cn_mean"]:
                i["freq"] = round(float(i["cn_mean"]), 2)
            elif "rnasv_reads" in i.keys() and i["rnasv_reads"]:
                i["freq"] = i["rnasv_reads"]
            else:
                i["freq"] = ""
            
            #融合拷贝数输出整数
            i['reads'] = '{}'.format(int(float(i['reads']))) if i.get('reads', '') else ''
            i['rnasv_reads'] = '{}'.format(int(float(i['rnasv_reads']))) if i.get('rnasv_reads', '') else ''

            '''
            # 用CARCINOGENICITY_GENE_LIST(致癌性)和FUNCTION_GENE_LIST(致病性)判定部分基因的解读等级，嵇梦晨，2024年2月1日
            if i['clinical_significance'] != '-' or i['function_classification'] != '-':
                if i['gene_symbol'] in CARCINOGENICITY_GENE_LIST:
                    i['clinic_g'] = clinic_strans_s.get(i['function_classification'], 3)
                    i['clinic_s'] = clinic_strans_s.get(i['function_classification'], 3)
                elif i['gene_symbol'] in FUNCTION_GENE_LIST:
                    i['clinic_g'] = clinic_strans_g.get(i['clinical_significance'], 3)
                    i['clinic_s'] = clinic_strans_g.get(i['clinical_significance'], 3)
                else:
                    i['clinic_g'] = clinic_strans_g.get(i['clinical_significance'], 3) if i['clinical_significance'] != '-' else clinic_strans_s.get(i['function_classification'], 3)
                    i['clinic_s'] = clinic_strans_s.get(i['function_classification'], 3) if i['function_classification'] != '-' else clinic_strans_g.get(i['clinical_significance'], 3)
            else:
                i['clinic_g'] = 3
                i['clinic_s'] = 3
            '''
            # snvindel返回两个解读，胚系看致病性体系看致癌性，返回一个就用对应的一个，CNV和SV不变
            if i["bio_category"] in ["Cnv", "Sv", "PSeqRnaSv"]:
                   i["clinic_g"] = clinic_strans_g.get(i["function_classification"])
                   i["clinic_s"] = clinic_strans_s.get(i["function_classification"])
                   i["clinic_num_g"] = clinic_strans_g_num.get(i["function_classification"], 3)
                   i["clinic_num_s"] = clinic_strans_s_num.get(i["function_classification"], 3)
            else:
                if i["clinical_significance"] != "-" or i["function_classification"] != "-":
                    i["clinic_g"] = clinic_strans_g.get(i["clinical_significance"]) if i["clinical_significance"] != "-" else clinic_strans_g.get(i["function_classification"])
                    i["clinic_s"] = clinic_strans_s.get(i["function_classification"]) if i["function_classification"] != "-" else clinic_strans_s.get(i["clinical_significance"])
                    i["clinic_num_g"] = clinic_strans_g_num.get(i["clinical_significance"], 3) if i["clinical_significance"] != "-" else clinic_strans_g_num.get(i["function_classification"], 3)
                    i["clinic_num_s"] = clinic_strans_s_num.get(i["function_classification"], 3) if i["function_classification"] != "-" else clinic_strans_s_num.get(i["clinical_significance"], 3)

            # CNV变异类型 可改为拷贝数变异，检测结果维持 转录本+扩增不变
            i["type_cn"] = "拷贝数变异" if i["bio_category"] == "Cnv" else "融合" if i["bio_category"] == "PSeqRnaSv" else "非移码突变" \
                if i["type"] == "nonSynonymous_Substitution" and re.search("del|ins", i["hgvs_p"]) else \
                type_strans.get(i["type"])
            
            # CNV增加转录本
            #i['cnv_transcript'] = cnv_transcript_cp40.get(i['gene_symbol']) if i['bio_category'] == 'Cnv' else ''
            #不确定升级后的流程是否输出CNV转录本, 2026年1月15日
            if i['bio_category'] == 'Cnv':
                if i["transcript_primary"] and i["transcript_primary"] != "-":
                    i['cnv_transcript'] = i["transcript_primary"]
                else:
                    i['cnv_transcript'] = cnv_transcript_cp40.get(i['gene_symbol'], '')

            # 处理药物，该字段放在用药提示部分
            i["regimen"] = []
            # 删除"KRAS", "HRAS", "NRAS", "PIK3CA", "FGFR4", "TP53"的诊断证据
            if 'evi_sum' in i.keys() and i['evi_sum']:
                #2026年1月14日 模板更新，保留TP53和PIK3CA的诊断证据
                if i['gene_symbol'] in ["KRAS", "HRAS", "NRAS", "FGFR4"]:
                #if i['gene_symbol'] in ["KRAS", "HRAS", "NRAS", "PIK3CA", "FGFR4", "TP53"]:
                    i['evi_sum'] = [item for item in i['evi_sum'] if item['evidence_type'] != 'Diagnostic']
            
            if "evi_sum" in i.keys() and i["evi_sum"]:
                for regimen_item in i["evi_sum"]:                   
                    if regimen_item["evidence_type"] == "Prognostic":
                        regimen_item["regimen_name"] = regimen_item["tumor_name_cn"] + "中预后较差" if regimen_item["clinical_significance"] == "Poor" else regimen_item["tumor_name_cn"] + "中预后较好"
                    if regimen_item["evidence_type"] == "Diagnostic":
                        if regimen_item["tumor_name_cn"]:
                            regimen_item["regimen_name"] = regimen_item["tumor_name_cn"] + "辅助诊断标志物"

                    regimen_item["sensi"] = sensi_stran.get(regimen_item["clinical_significance"], "0")
                    regimen_item["level"] = regimen_item["evi_conclusion"][0] if regimen_item["evi_conclusion"] else "0"
                    regimen_item["regimen_name_py"] = self._topinyin(regimen_item["regimen_name"]) if regimen_item["regimen_name"] else "0"
                    regimen_item["sense_rule"] = self.sense_rule(regimen_item["clinical_significance"])
                    regimen_item["regimen_summary"] = regimen_item["regimen_name"] + "(" + regimen_item["sensi"] + "，" + \
                                                      regimen_item["level"] + "级" + ")"
                    regimen_item = self.get_appr_info(regimen_item, regimen_dict, regimen_adaptation)

                    if "NMPA" in regimen_item["evi_origin_ZJZL"]:
                        regimen_item["judge_NMPA"] = 0
                    else:
                        regimen_item["judge_NMPA"] = 1
                    # 添加特殊药物的index，便于后续EGFR药物排序
                    regimen_item["regimen_index"] = regimen_rule.index(regimen_item["regimen_name"]) if regimen_item["regimen_name"] in regimen_rule else len(regimen_rule)
                
                # 和CP一致 2024年3月11日
                i["evi_sum"] = sorted(i["evi_sum"], key=lambda j : (j["level"], j["sense_rule"], j["judge_NMPA"], j["regimen_name_py"].upper()))
                if i['gene_symbol'] == 'EGFR':
                    i['evi_sum'] = sorted(i["evi_sum"], key=lambda j : j["regimen_index"])
            # 不存在用药
            else:
                pass

            # regimen存储用药信息
            i["regimen"] = [
                {
                    'regimen_name' : j['regimen_name'],
                    'appr_note' : j['appr_note'],
                    'sensi' : j['sensi'],
                    'level' : j['level'],
                    'evi_origin_ZJZL' : j['evi_origin_ZJZL'],
                    'evi_interpretation' : j['evi_interpretation']
                }
                for j in i["evi_sum"] if j["regimen_name"] and not re.search("预后", j["regimen_name"]) and not re.search("诊断", j["regimen_name"])
            ]
            i["regimen"] = i["regimen"] if i["regimen"] else []
                                                               
            i["Prognostic"] = [
                {
                    'tumor_name_cn' : j['tumor_name_cn'],
                    'evi_interpretation' : j['evi_interpretation']
                }
                for j in i["evi_sum"] if j["regimen_name"] and re.search("预后", j["regimen_name"])
            ]

            i["Diagnostic"] = [
                {
                    'tumor_name_cn' : j['tumor_name_cn'],
                    'evi_interpretation' : j['evi_interpretation']
                }
                for j in i["evi_sum"] if j["regimen_name"] and re.search("诊断", j["regimen_name"])
            ]
            
            # 去重
            i['Diagnostic'] = self.merge_diag(i['Diagnostic'])
            i['Prognostic'] = self.merge_diag(i['Prognostic'])

            i["Prognostic_str"] = [j["regimen_name"] for j in i["evi_sum"] if j["regimen_name"] and re.search("预后", j["regimen_name"])]
            i["Diagnostic_str"] = [j["regimen_name"] for j in i["evi_sum"] if j["regimen_name"] and re.search("诊断", j["regimen_name"])]
            i["Prognostic_str"] = list(set(i["Prognostic_str"]))
            i["Prognostic_str"] = '\n'.join(i["Prognostic_str"] if i["Prognostic_str"] else ["-"])
            i["Diagnostic_str"] = list(set(i["Diagnostic_str"]))
            i["Diagnostic_str"] = '\n'.join(i["Diagnostic_str"] if i["Diagnostic_str"] else ["-"]) 

            # 1. 添加体细胞/胚系标签，便于排序 体细胞 > 胚系
            i["var_ori_level"] = 2 if i["var_origin"] == "germline" else 1

            # 汇总用药等级，判断为明确或是潜在变异
            # 不含遗传易感的证据
            i['evi_sum'] = [k for k in i['evi_sum'] if k['evidence_type'] != 'Predisposing'] if "evi_sum" in i.keys() and i["evi_sum"] else []

            level_list = [regimen_item["level"] for regimen_item in i["evi_sum"]] if "evi_sum" in i.keys() and i["evi_sum"] else []
            i['top_level'] = 4 if "A" in level_list else 3 if "B" in level_list else 2 if "C" in level_list else 1 if "D" in level_list else 0
            
            #新增-变异类型合集，包含用药、诊断和预后三类 2024年12月19日
            evi_type_list = set([item["evidence_type"] for item in i["evi_sum"]])
            i["evi_type_list"] = [i for i in evi_type_list if i in ["Predictive", "Prognostic", "Diagnostic"]]
            evitype_dict = {"Predictive" : "用药", "Prognostic" : "预后", "Diagnostic" : "诊断"}
            i["evi_type_list_str"] = "、".join([evitype_dict.get(i, i) for i in i["evi_type_list"]])+"相关" if i["evi_type_list"] else ""

            # level_str用于变异排序
            i["level_str"] = "I类" if set(level_list) & {"A", "B"} else "II类" if set(level_list) & {"C", "D"} \
                else "II类" if i["clinic_g"] in ["致病性变异", "疑似致病性变异"] and i["var_origin"] == "germline" else "III类"
            

            # 添加变异类型标签，便于排序 snvindel > sv > cnv
            i['var_type_num'] = 3 if i['bio_category'] == 'Snvindel' else 2 if i['bio_category'] == 'PSeqRnaSv' else 1 if i['bio_category'] == 'Cnv' else 0
            
            # 避免基因功能输出None,jmc, 2024年9月6日
            if i['bio_category'] == 'PSeqRnaSv':
                i['three_prime_gene_function'] = '' if not i['three_prime_gene_function'] else i['three_prime_gene_function']
                i['five_prime_gene_function'] = '' if not i['five_prime_gene_function'] else i['five_prime_gene_function']

            var_list.append(copy.deepcopy(i))

        # 特殊变异：RNA MET融合和DNA MET 14跳跃
        judge_MET14 = []
        judge_METSV = []
        for var_item in var_list:
            if var_item["bio_category"] == "PSeqRnaSv" and var_item["five_prime_gene"] == "MET" and \
                    var_item["three_prime_gene"] == "MET":
                judge_METSV = var_item
                regimen_METSV = var_item["regimen"]
                var_item["MET14"] = "yes"
        if judge_METSV:
            for var_item in var_list:
                if var_item["bio_category"] == "Snvindel" and var_item["gene_symbol"] == "MET" and \
                        regimen_METSV == var_item["regimen"]:
                    judge_MET14 = var_item
                    var_list.remove(var_item)
                    var_list.remove(judge_METSV)
                    judge_MET14["hgvs_p_2"] = "MET-MET融合"
                    judge_MET14["freq2"] = str(int(judge_METSV["freq"])) + "copies"
        if judge_MET14:
            var_list.insert(0, judge_MET14)
        # 特殊变异结束

        # 融合变异特殊处理
        # 保证报出的基因为1个
        funcgene_list = ['ALK', 'BCOR', 'BRAF', 'CIC', 'EWSR1', 'FOSB', 'FUS', 'GLI1', 'MIR143', 'NR4A3',  'PAX3',  'PHF1', 'RAF1', 'RET', 'ROS1', 'SS18', 'USP6', 'YAP1', 'YWHAE',
        'CSF1', 'ESR1', 'GREB1', 'NCOA2', 'NCOA3', 'NTRK1', 'NTRK2', 'NTRK3', 'PGR', 'PLAG1', 'RAD51B', 'TFE3']
        for var in var_list:
            if var['bio_category'] == 'PSeqRnaSv' and ',' in var['gene_symbol']:
                if var['gene_symbol'] == 'EWSR1,NR4A3':
                    var['gene_symbol'] = 'EWSR1'
                else:
                    five_prime_gene, three_prime_gene = var['gene_symbol'].split(',')
                    if not (five_prime_gene in funcgene_list and three_prime_gene in funcgene_list):
                        var['gene_symbol'] = three_prime_gene
                    elif five_prime_gene in funcgene_list and three_prime_gene in funcgene_list:
                        var['gene_symbol'] = three_prime_gene
                    elif three_prime_gene in funcgene_list:
                        var['gene_symbol'] = three_prime_gene
                    else:
                        var['gene_symbol'] = five_prime_gene
                    
        # var排序
        var_list = sorted(var_list, key=lambda j: (
            j["var_ori_level"], j['level_str'], j['top_level'], j['clinic_s'], j["var_type_num"], float(str(j["freq"]).strip("%"))), reverse=True)
        
        return var_list
    
    # 获取治疗方案的证据来源/获批机构，和适应症 2024年3月7日
    def get_appr_info(self, evi, regimen_dict, regimen_adaptation):
        evi["regimen_summary"] = "" 
        evi["appr_info"] = [] 
	    
        apprlist = ["FDA", "NMPA", "NCCN", "CSCO"]
        appdict = {"FDA" : 0, "NCCN" : 1, "NMPA" : 2, "CSCO" : 3}

        if re.search('Sensitive', evi['clinical_significance']):
            regimen_refer_agency_ZJZL = regimen_dict.get(evi["regimen_name"], "")
        else:
            regimen_refer_agency_ZJZL = '/'.join(sorted(list(set(apprlist) & set(re.split(",", evi["refer_agency"]))),key=lambda i:appdict.get(i))) if "refer_agency" in evi.keys() and evi["refer_agency"] and set(apprlist) & set(re.split(",", evi["refer_agency"])) else ""
        
        evidence_level = evi_level_dict.get(evi["evidence_level"], evi["evidence_level"]) if "evidence_level" in evi.keys() and evi["evidence_level"] else ""
        
        evi['evi_origin_ZJZL'] = ''
        evi['appr_note'] = ''
        if evi['level'] == 'A':
            evi["evi_origin_ZJZL"] = regimen_refer_agency_ZJZL if regimen_refer_agency_ZJZL else evidence_level if evidence_level else "-"
            if regimen_refer_agency_ZJZL and evi["regimen_name"] and evi["regimen_name"] in regimen_adaptation.keys() and regimen_adaptation[evi["regimen_name"]]["adaptation"]:
                evi["appr_note"] = "yes"
                evi["appr_info"].append({
                    'regimen_cn' : regimen_adaptation[evi["regimen_name"]]["regimen_cn"],
                    'regimen_en' : regimen_adaptation[evi["regimen_name"]]["regimen_en"],
                    'appr_info' : regimen_adaptation[evi["regimen_name"]]["adaptation"]
                })
        elif evi['level'] == 'C':
            if evi["evi_conclusion"] == "C3":
                if evidence_level and evidence_level == "临床试验":
                    evi["evi_origin_ZJZL"] = "其他癌种中获批/临床试验"
                else:
                    evi["evi_origin_ZJZL"] = "其他癌种中获批"
            else:
                evi["evi_origin_ZJZL"] = evidence_level if evidence_level else "-"
            
            if evi["evi_conclusion"] == "C3" and regimen_refer_agency_ZJZL and evi["regimen_name"] and  evi["regimen_name"] in regimen_adaptation.keys() and regimen_adaptation[evi["regimen_name"]]["adaptation"]:
                evi["appr_note"] = "yes"
                evi["appr_info"].append({
                    'regimen_cn' : regimen_adaptation[evi["regimen_name"]]["regimen_cn"],
                    'regimen_en' : regimen_adaptation[evi["regimen_name"]]["regimen_en"],
                    'appr_info' : regimen_adaptation[evi["regimen_name"]]["adaptation"]
                })
        elif evi['level'] in ['B', 'D']:
            evi["evi_origin_ZJZL"] = evidence_level if evidence_level else "-"
        else:
            evi["evi_origin_ZJZL"] = "-"
        
        if evi['appr_note']:
            evi["regimen_summary"] = evi["regimen_name"]+"※（"+evi["sensi"]+"，"+evi["level"]+"级，"+evi["evi_origin_ZJZL"]+"）"
        else:
            evi["regimen_summary"] = evi["regimen_name"]+"（"+evi["sensi"]+"，"+evi["level"]+"级，"+evi["evi_origin_ZJZL"]+"）"
        
        return evi
    
    def get_regimen_appr(self, regimen_info):
        regimen_dict = {}
        regimen_adaptation = {}
        if regimen_info:
            for regimen in regimen_info:
                apprlist = ["FDA", "NMPA", "NCCN", "CSCO"]
                appdict = {"FDA" : 0, "NCCN" : 1, "NMPA" : 2, "CSCO" : 3}
                regimen_appr = '/'.join(sorted(list(set(apprlist) & set(regimen["approval_organization"])),key=lambda i:appdict.get(i))) if regimen["approval_organization"] and set(apprlist) & set(regimen["approval_organization"]) else '-'
                regimen_dict[regimen["regimen_cn"]] = regimen_appr
                regimen_dict[regimen["regimen_en"]] = regimen_appr

                adaptation = list(itertools.chain(*[re.split("\n", i.strip()) for i in regimen["adaptation_disease_cn"]])) if regimen["adaptation_disease_cn"] else []
                if adaptation:
                    regimen_adaptation[regimen["regimen_cn"]] = {
                        "regimen_cn" : regimen["regimen_cn"],
                        "regimen_en" : regimen["regimen_en"],
                        "adaptation" : "".join(adaptation)
                    }
                    regimen_adaptation[regimen["regimen_en"]] = {
                        "regimen_cn" : regimen["regimen_cn"],
                        "regimen_en" : regimen["regimen_en"],
                        "adaptation" : "".join(adaptation)
                    }
        return regimen_dict, regimen_adaptation


    def _merge_evi(self, var_list):
        # 定义合并规则： 字段名，连接符，取值函数
        MERGE_RULES = [
            ('regimen_name_list', '、', lambda x: x["regimen_name"]),
            ('evi_conclusion_list', '/', lambda x: x["evi_conclusion"][0] if x["evi_conclusion"] else ""),
            ('tag_list', '/', lambda x: str(x['tag'])),
            ('regimen_py_list', '/', lambda x: x["regimen_name"]),
        ]
        for var in var_list:
            if not var.get("evi_sum"):
                var['evi_sum_merge'] = []
                continue

            groups = defaultdict(list)
            for item in var["evi_sum"]:
                groups[item["evi_interpretation"]].append({
                    'regimen_name': item["regimen_name"],
                    'evi_conclusion': item["evi_conclusion"],
                    'tag': sensi_stran.get(item["clinical_significance"]),
                    'regimen_name_py': self._topinyin(item["regimen_name"]),
                })
            # 构建合并结果
            merged = []
            for interpretation, items in groups.items():
                item_dict = {"evi_interpretation": interpretation}
                for field, delimiter, get_value in MERGE_RULES:
                    values = [get_value(item) for item in items]
                    # 特殊处理空值过滤
                    if field == "evi_conclusion_list":
                        values = [v for v in values if v]
                    item_dict[field] = delimiter.join(filter(None, values))
                merged.append(item_dict)
            # 排序逻辑
            var["evi_sum_merge"] = sorted(
                merged,
                key=lambda x: (
                    x["evi_conclusion_list"],
                    x["tag_list"],
                    x["regimen_py_list"]
                )
            )

            # 优化类型判断逻辑
            evidence_types = {i["evidence_type"] for i in var["evi_sum"]}
            var["evi_type"] = (
                "Diagnostic" if evidence_types == {"Diagnostic"} else
                "Prognostic" if evidence_types == {"Prognostic"} else
                "Predictive"
            )

        return var_list
    
    def var_inter(self, var_list):
        evitype_dict = {"Predictive" : "用药相关", "Prognostic" : "预后相关", "Diagnostic" : "诊断相关"}
        level_dict = {"A" : "临床意义明确", "B" : "临床意义明确", "C" : "有潜在临床意义", "D" : "有潜在临床意义"}
        level_priority = {'A': 4, 'B': 3, 'C': 2, 'D': 1}
        
        for var in var_list:
            type_max_level = {}
            for evi in var.get('evi_sum', []):
                evi_type = evi.get('evidence_type', '')
                level = evi.get('level', '')

                if not evi_type or level not in level_dict.keys():
                    continue

                current_priority = level_priority[level]
                if evi_type not in type_max_level.keys() or current_priority > level_priority[type_max_level[evi_type]]:
                    type_max_level[evi_type] = level
            
            grouped_results = {
                '临床意义明确': set(),
                '有潜在临床意义': set()
            }

            for evi_type, max_level in type_max_level.items():
                group = level_dict[max_level]
                translated_type = evitype_dict.get(evi_type, evi_type)
                grouped_results[group].add(translated_type)
            
            output = []
            for group in ['临床意义明确', '有潜在临床意义']:
                types = grouped_results[group]
                if types:
                #if types := grouped_results[group]:
                    sorted_types = '/'.join(sorted(types))
                    output.append(f"{group}（{sorted_types}）")
            
            var['zjzl_type_str'] = '\n'.join(output) if output else '临床意义不明确'

        return var_list
        
    def _var_count(self):
        """
        统计变异个数
        :return:
        """
        var = self._var_list()

        somatic_all = len([i for i in var if self.origin(i) == "S"])
        # 体系
        somatic_I = len([i for i in var if self.origin(i) == "S" and i["clinic_s"] == 5])
        somatic_II = len([i for i in var if self.origin(i) == "S" and i["clinic_s"] == 4])
        somatic_III = len([i for i in var if i["clinic_s"] == 3])
        somatic_drug = len([i for i in var if self.origin(i) == "S" and 'evi_sum' in i.keys() and i['evi_sum']])
        
        diagnostic = len([i for i in var if i['Diagnostic']])
        
        # 区分靶向用药、诊断、预后计数
        predictive = len([i for i in var if i['regimen']])

        # 胚系
        germline_4_5 = len([i for i in var if self.origin(i) == "G" and i["clinic_g"] in [5, 4]])
        germline_drug = len([i for i in var if self.origin(i) == "G" and "evi_sum" in i.keys() and i["evi_sum"]])
        return {'somatic_all': somatic_all, 'somatic_I': somatic_I, 'somatic_II': somatic_II,
                'somatic_III': somatic_III, 'somatic_drug' : somatic_drug, 
                'germline_4_5': germline_4_5, 'germline_drug': germline_drug, 'diagnostic': diagnostic, 'predictive' : predictive}


    def get_base_info(self):
        """
        个人基本信息
        :return:
        """
        # 样本id处理
        self.data_js['result']['sample_info']['sample_id'] = self._sample_id()
        # 年龄
        age = self.data_js['sample_info'].get('age', '')
        self.data_js['result']['sample_info']['age'] = str(age) + ' 岁' if age else '-'
        # 报告日期（当日）
        self.data_js['result']['sample_info']['report_date'] = datetime.datetime.now().strftime("%Y年%m月%d日")

        # 分子特检号
        molecular_number = self.data_js['sample_info'].get('molecular_number', '-')
        self.data_js['result']['sample_info']['molecular_number'] = molecular_number

        return self.data_js['result']['sample_info']

    def get_var(self):
        """
        变异结果汇总
        :return: self.data_js['result']['var']['level_I']
        :return: self.data_js['result']['var']['level_II']
        :return: self.data_js['result']['var']['level_III']
        """
        var_li = self._var_list()
        var_li = self._merge_evi(var_li)
        var_li = self.var_inter(var_li)
        var_li = sorted(var_li, key=lambda j: (
            j["var_ori_level"], j['level_str'], j['top_level'], j['clinic_s'], j["var_type_num"], float(str(j["freq"]).strip("%"))), reverse=True)

        self.data_js['result']['var']['all_var'] = var_li
        self.data_js['result']['var']['level_I'] = [i for i in var_li if i['level_str'] == 'I类']
        self.data_js['result']['var']['level_II'] = [i for i in var_li if i['level_str'] == 'II类']
        self.data_js['result']['var']['level_III'] = [i for i in var_li if i['level_str'] == 'III类']
        return self.data_js['result']['var']
    
    def get_test_result(self):
        """
        结果小结
        :return: self.data_js['result']['var_count']
        """
        var_num_dict = self._var_count()

        self.data_js['result']['var_count'] = var_num_dict
        return self.data_js['result']['var_count']
    
    SARCOMA_RULES = {
        "H3_result" : {
            "type" : "snvindel",
            "genes" : ["H3-3A"]
        },
        "CDK4_result" : {
            "type" : "cnv",
            "genes" : ["CDK4", "MDM2"]
        },
        "DDIT3_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("FUS", "DDIT3"), ("EWSR1", "DDIT3")]
        },
        "NAB2_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("NAB2", "STAT6")]
        },
        "ALK_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("TPM3", "ALK"), ("TPM4", "ALK"), ("CLTC", "ALK"), ("RANBP2", "ALK"),
                ("ATIC", "ALK"), ("CARS1", "ALK"), ("SEC31A", "ALK"), ("PPFIBP1", "ALK"),
                ("TIMP3", "ALK"), ("IGFBP5", "ALK"), ("THBS1", "ALK"), ("RRBP1", "ALK"),
                ("ETV6", "NTRK3"), ("TFG", "ROS1"), ("YWHAE", "ROS1")
            ]
        },
        "COL1A1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("COL1A1", "PDGFB"), ("COL6A3", "PDGFD"), ("EMILIN2", "PDGFD")]
        },
        "ETV6_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("ETV6", "NTRK3")]
        },
        "FUS_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("FUS", "CREB3L1"), ("FUS", "CREB3L2")]
        },
        "EWSR1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("EWSR1", "CREB3L1"), ("FUS", "CREB3L2"),
                ("FUS", "CREB3L1"), ("YAP1", "KMT2A")
            ]
        },
        "ATF1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("EWSR1", "ATF1"), ("FUS", "ATF1"), ("EWSR1", "CREB1")]
        },
        "CTNNB1_result" : {
            "type" : "snvindel",
            "genes" : ["CTNNB1"]
        },
        "CSF1_result" : {
            "type" : "fusion_any",
            "genes" : ["CSF1"]
        },
        "WWTR1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("WWTR1", "CAMTA1"), ("YAP1", "TFE3"), ("MBNL1", "FOS"),
                ("VIM", "FOS"), ("ZFP36", "FOSB"), ("WWTR1", "FOSB"),
                ("ACTB", "FOSB"), ("SERPINE1", "FOSB"), ("YAP1", "MAML2"),
                ("PTBP1", "MAML2")
            ]
        },
        "MYC_result" : {
            "type" : "cnv",
            "genes" : ["MYC"]
        },
        "DICER1_result" : {
            "type" : "snvindel",
            "genes" : ["DICER1"]
        },
        "PAX3_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("PAX3", "FOXO1"), ("PAX7", "FOXO1"), ("PAX3", "FOXO4"),
                ("PAX3", "NCOA1"), ("PAX3", "NCOA2"), ("FOXO1", "FGFR1"),
                ("PAX3", "INO80D")
            ]
        },
        "RMS_result" : {
            "type" : "snvindel",
            "genes" : ["PIK3CA", "TP53"]
        },
        "NCOA2_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("SRF", "NCOA2"), ("TEAD1", "NCOA2"), ("VGLL2", "NCOA2"),
                ("VGLL2", "CITED2")
            ]
        },
        "HEY1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("HEY1", "NCOA2")]
        },
        "NF1_result" : {
            "type" : "snvindel",
            "genes" : ["NF1", "CDKN2A"]
        },
        "ZNF444_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("EWSR1", "POU5F1"), ("EWSR1", "PBX1"),
                ("FUS", "KLF17"), ("EWSR1", "PBX3"),
                ("EWSR1", "ZNF444")
            ]
        },
        "NTRK_result" : {
            "type" : "fusion_any",
            "genes" : ["NTRK1", "NTRK2", "NTRK3"]
        },
        "SS18_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("SS18", "SSX1"), ("SS18", "SSX2"),
                ("SS18", "SSX4"), ("SS18L1", "SSX1")
            ]
        },
        "SMARCB1_result" : {
            "type" : "snvindel",
            "genes" : ["SMARCB1"]
        },
        "ASPSCR1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("ASPSCR1", "TFE3")]
        },
        "CREB1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("EWSR1", "ATF1"), ("EWSR1", "CREB1")]
        },
        "IDH1_result" : {
            "type" : "snvindel",
            "genes" : ["IDH1", "IDH2"]
        },
        "NR4A3_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("EWSR1", "NR4A3"), ("TAF15", "NR4A3"),
                ("TCF12", "NR4A3"), ("TFG", "NR4A3")
            ]
        },
        "WT1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("EWSR1", "WT1")]
        },
        "BCOR_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("CIC", "DUX4"), ("BCOR", "CCNB3")
            ]
        },
        "TFE3_result1" : {
            "type" : "fusion_any",
            "genes" : ["TFE3", "RAD51B"]
        },
        "TFE3_result2" : {
            "type" : "snvindel",
            "genes" : ["TSC1", "TSC2"]
        },
        "TFE3_result3" : {
            "type" : "fusion_ordered",
            "pairs" : [("HTR4", "ST3GAL1")]
        },
        "APC_result" : {
            "type" : "snvindel",
            "genes" : ["APC", "CTNNB1"]
        },
        "ETV1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("EWSR1", "FLI1"), ("EWSR1", "ERG"),
                ("EWSR1", "FEV"), ("EWSR1", "ETV1"),
                ("EWSR1", "ETV4"), ("EWSR1", "PATZ1"),
                ("FUS", "ERG"), ("FUS", "FEV")
            ]
        },
        "CIC_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("CIC", "DUX4"), ("CIC", "FOXO4"),
                ("CIC", "NUTM1"), ("CIC", "NUTM2B")
            ]
        },
        "CCNB3_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("BCOR", "CCNB3"), ("BCOR", "MAML3"),
                ("ZC3H7B", "BCOR"), ("YWHAE", "NUTM2B")
            ]
        },
        "YWHAE_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("YWHAE", "NUTM2B")]
        },
        "SP3_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("EWSR1", "NFATC2"), ("EWSR1", "SP3"),
                ("EWSR1", "POU5F1"), ("EWSR1", "PATZ1"),
                ("EWSR1", "SMARCA5"), ("FUS", "NFATC2")
            ]
        },
        "SQSTM1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("SQSTM1", "ALK"), ("VCL", "ALK")]
        },
        "TFCP2_result" : {
            "type" : "fusion_ordered",
            "pairs" : [
                ("EWSR1", "TFCP2"), ("FUS", "TFCP2"),
                ("MEIS1", "NCOA2")
            ]
        },
        "GLI1_result1" : {
            "type" : "fusion_ordered",
            "pairs" : [("ACTB", "GLI1"), ("MALAT1", "GLI1"), ("PTCH1", "GLI1")]
        },
        "GLI1_result2" : {
            "type" : "cnv",
            "genes" : ["GLI1"]
        },
        "MALAT1_result" : {
            "type" : "fusion_ordered", 
            "pairs" : [("MALAT1", "GLI1")]
        },
        "SRF_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("SRF", "ICA1L")]
        },
        "OGT_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("OGT", "FOXO3")]
        },
        "MEIS1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("MEIS1", "NCOA2")]
        },
        "CREB2_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("EWSR1", "CREB1")]
        },
        "MIR143_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("MIR143", "NOTCH1"), ("MIR143", "NOTCH2"), ("MIR143", "NOTCH3")]
        },
        "RAF1_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("PDZRN3", "RAF1"), ("SLMAP", "RAF1"), ("TMF1", "RAF1"), ("MTAP", "RAF1"), ("TFG", "RET"), ("MYH10", "RET"),
            ("NCOA4", "RET"), ("VCL", "RET"), ("CLIP2", "RET"), ("CCDC6", "RET"), ("KHDRBS1", "RET"), ("SPECC1L", "RET"),
            ("KIAA1217", "RET"), ("PPP1CB", "ALK"), ("CUX1", "BRAF"), ("SEPTIN7", "BRAF"), ("CDC42SE2", "BRAF")]
        },
        "PGR_result" : {
            "type" : "fusion_any",
            "genes" : ["PGR"]
        },
        "PLAG1_result" : {
            "type" : "fusion_any",
            "genes" : ["PLAG1"]
        },
        "USP6_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("CDH11", "USP6"), ("CNBP", "USP6"), ("COL1A1", "USP6"), ("MYH9", "USP6"),
            ("OMD", "USP6"), ("THRAP3", "USP6")]
        },
        "NUTM2A_result" : {
            "type" : "fusion_ordered",
            "pairs" : [("YWHAE", "NUTM2A"), ("YWHAE", "NUTM2B"), ("ZC3H7B", "BCOR")]
        },
        "JAZF1_result" : {
            "type" : "fusion_ordered",
            "pairs": [("JAZF1","SUZ12"), ("JAZF1","PHF1"), ("MEAF6","PHF1"),
            ("EPC1","PHF1"), ("MBTD1","EZHIP"), ("BRD8","PHF1"), ("EPC2","PHF1"), ("EPC1","SUZ12")]
        },
        "ESR1_result" : {
            "type" : "fusion_any",
            "genes" : ["ESR1", "GREB1"]
        },
        "NCOA3_result" : {
            "type" : "fusion_any",
            "genes" : ["NCOA3", "NCOA2"]
        },
        "SMARCA4_result" : {
            "type" : "snvindel",
            "genes" : ["SMARCA4"]
        },
        "SDH_result" : {
            "type" : "fusion_any",
            "genes" : ["SDHA", "SDHB", "SDHC", "SDHD"]
        },
        "NUTM1_result" : {
            "type" : "fusion_any",
            "genes" : ["NUTM1"]
        }
    }

    def _match_sarcoma_rule(self, var, rule):
        bio = var.get("bio_category")
        if rule["type"] == "fusion_ordered":
            if bio != "PSeqRnaSv":
                return False
            five = var.get("five_prime_gene")
            three = var.get("three_prime_gene")
            return (five, three) in rule.get("pairs", [])
        if rule["type"] == "fusion_any":
            if bio != "PSeqRnaSv":
                return False
            five = var.get("five_prime_gene")
            three = var.get("three_prime_gene")
            genes = rule.get("genes", set())
            return five in genes or three in genes
        if rule["type"] == "cnv":
            return (
                bio == "Cnv"
                and var.get("gene_symbol") in rule.get("genes", set())
            )
        if rule["type"] == "snvindel":
            return (
                bio == "Snvindel"
                and var.get("gene_symbol") in rule.get("genes", set())
            )
        return False

    def dedup(self, li):
        new = []
        for var in li:
            if var not in new:
                new.append(var)
        return new
    
    def get_sarcoma_result(self):
        results = {k: [] for k in self.SARCOMA_RULES.keys()}
        for var in self._var_list():
            for result_key, rule in self.SARCOMA_RULES.items():
                if self._match_sarcoma_rule(var, rule):
                    if var["bio_category"] == "PSeqRnaSv":
                        results[result_key].append(
                            f'{var["five_prime_gene"]}-{var["three_prime_gene"]}'
                        )
                    elif var["bio_category"] == "Cnv":
                        results[result_key].append(var["gene_symbol"])
                    else:
                        results[result_key].append(var.get("gene_symbol"))
        for k, v in results.items():
            v = self.dedup(v)
            self.data_js["result"][k] = list(dict.fromkeys(v))
        
        return
    
    def get_rouliu_result(self):
        """
        肉瘤指南或共识推荐标志物检测结果汇总
        :return:
        """
        NTRK_fusion = ["NTRK1", "NTRK2", "NTRK3"]
        CDK4_cnv = ["CDK4"]
        MYC_cnv = ["MYC"]
        DDIT3_fusion = ["FUS", "DDIT3", "EWSR1"]
        NAB2_fusion = ["NAB2", "STAT6"]
        ALK_fusion = ["TPM3", "ALK", "TPM4", "CLTC", "RANBP2", "ATIC", "CARS", "SEC31L1", "PPFIBP1", "TIMP3", "IGFBP5", "THBS1", "RRBP1", "ETV6", "NTRK3", "TFG", "ROS1", "YWHAE"]
        COL1A1_fusion = ["COL1A1", "PDGFB", "COL6A3", "PDGFD", "EMILIN2"]
        ETV6_fusion = ["ETV6", "NTRK3"]
        FUS_fusion = ["FUS", "CREB3L1", "CREB3L2"]
        EWSR1_fusion = ["EWSR1", "CREB3L1", "FUS", "CREB3L2", "YAP1", "KMT2A"]
        ATF1_fusion = ["EWSR1", "ATF1", "CREB1", "FUS"]
        CSF1_fusion = ["CSF1"]
        WWTR1_fusion = ["WWTR1", "CAMTA1", "YAP1", "TFE3", "MBNL1", "FOS", "VIM", "ZFP36", "FOSB", "ACTB", "SERPINE1", "MAML2", "PTBP1", "MAML2"]
        PAX3_fusion = ["PAX3", "FOXO1", "PAX7", "FOXO4", "NCOA1", "NCOA2", "FGFR1", "INO80D"]
        NCOA2_fusion = ["SRF", "NCOA2", "TEAD1", "NCOA2", "VGLL2", "NCOA2", "VGLL2", "CITED2"]
        HEY1_fusion = ["HEY1", "NCOA2"]
        ZNF444_fusion = ["EWSR1", "POU5F1", "PBX1", "FUS", "KLF17", "PBX3", "EWSR1", "ZNF444"]
        SS18_fusion = ["SS18", "SSX1", "SSX2", "SSX4", "SS18L1"]
        ASPSCR1_fusion = ["ASPSCR1", "TFE3"]
        CREB1_fusion = ["EWSR1", "ATF1", "CREB1"]
        NR4A3_fusion = ["EWSR1", "NR4A3", "TAF15", "TCF12", "TFG"]
        WT1_fusion = ["EWSR1", "WT1"]
        BCOR_fusion = ["CIC", "DUX4", "BCOR", "CCNB3"]
        TFE3_fusion = ["TFE3", "RAD51B", "HTR4", "ST3GAL1"]
        ETV1_fusion = ["EWSR1", "FLI1", "ERG", "FEV", "ETV1", "ETV4", "PATZ1", "FUS"]
        CIC_fusion = ["CIC", "DUX4", "FOXO4", "NUTM1", "NUTM2B"]
        CCNB3_fusion = ["BCOR", "CCNB3", "MAML3", "ZC3H7B", "YWHAE", "NUTM2B"]
        YWHAE_fusion = ["YWHAE", "NUTM2B"]
        SP3_fusion = ["EWSR1", "NFATC2", "SP3", "POU5F1", "PATZ1", "SMARCA5", "FUS", "NFATC2"]
        SQSTM1_fusion = ["SQSTM1", "ALK", "VCL"]
        TFCP2_fusion = ["EWSR1", "TFCP2", "FUS", "MEIS1", "NCOA2"]
        GLI1_fusion = ["ACTB", "GLI1", "MALAT1", "PTCH1"]
        MALAT1_fusion = ["MALAT1", "GLI1"]
        SRF_fusion = ["SRF", "ICA1L"]
        OGT_fusion = ["OGT", "FOXO3"]
        MEIS1_fusion = ["MEIS1", "NCOA2"]
        CREB2_fusion = ["EWSR1", "CREB1"]
        MIR143_fusion = ["MIR143", "NOTCH1", "NOTCH2", "NOTCH3"]
        RAF1_fusion = ["PDZRN3", "RAF1", "SLMAP", "TMF1", "MTAP", "TFG", "RET", "MYH10", "NCOA4", "VCL", "CLIP2", "CCDC6", "KHDRBS1", "SPECC1L", "KIAA1217", "PPP1CB", "ALK", "CUX1", "BRAF", "SEPTIN7", "CDC42SE2"]
        PGR_fusion = ["PGR"]
        PLAG1_fusion = ["PLAG1"]
        USP6_fusion = ["CDH11", "USP6", "CNBP", "COL1A1", "MYH9", "OMD", "THRAP3"]
        NUTM2A_fusion = ["YWHAE", "NUTM2A", "NUTM2B", "ZC3H7B", "BCOR"]
        JAZF1_fusion = ["JAZF1", "SUZ12", "PHF1", "MEAF6", "EPC1", "MBTD1", "CXorf67", "BRD8", "PHF1", "EPC2", "EPC1", "SUZ12"]
        ESR1_fusion = ["ESR1", "GREB1"]
        NCOA3_fusion = ["NCOA2", "NCOA3"]
        CTNNB1_gene  = ["CTNNB1"]
        #KRAS_gene = ["KRAS", "HRAS", "NRAS", "PIK3CA", "FGFR4"]
        #TP53_gene = ["TP53"]
        IDH1_gene = ["IDH1", "IDH2"]
    
        NTRK_result = []
        CDK4_result = []
        DDIT3_result = []
        NAB2_result = []
        ALK_result = []
        COL1A1_result = []
        ETV6_result = []
        FUS_result = []
        EWSR1_result = []
        ATF1_result = []
        CSF1_result = []
        WWTR1_result = []
        MYC_result = []
        PAX3_result = []
        NCOA2_result = []
        HEY1_result = []
        ZNF444_result = []
        SS18_result = []
        ASPSCR1_result = []
        CREB1_result = []
        NR4A3_result = []
        WT1_result = []
        BCOR_result = []
        TFE3_result = []
        ETV1_result = []
        CIC_result = []
        CCNB3_result = []
        YWHAE_result = []
        SP3_result = []
        SQSTM1_result = []
        TFCP2_result = []
        GLI1_result = []
        MALAT1_result = []
        SRF_result = []
        OGT_result = []
        MEIS1_result = []
        CREB2_result = []
        MIR143_result = []
        RAF1_result = []
        PGR_result = []
        PLAG1_result = []
        USP6_result = []
        NUTM2A_result = []
        JAZF1_result = []
        ESR1_result = []
        NCOA3_result = []
        CTNNB1_result = []
        #KRAS_result = []
        #TP53_result = []
        IDH1_result = []
        for i in self._var_list():
            var_info = "扩增" if i["bio_category"] == "Cnv" else i["five_prime_gene"] + "-" + i["three_prime_gene"] + "融合" if \
                i["bio_category"] == "PSeqRnaSv" else "" if i["bio_category"] == "knb" else i["hgvs_p"] if i["hgvs_p"] != "p.?" else i["hgvs_c"]
            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in NTRK_fusion or i["five_prime_gene"] in NTRK_fusion):
                NTRK_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])
                
            if i["bio_category"] == "Cnv" and i["gene_symbol"] in CDK4_cnv:
                CDK4_result = ["CDK4"]

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in DDIT3_fusion and i["five_prime_gene"] in DDIT3_fusion):
                DDIT3_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in NAB2_fusion and i["five_prime_gene"] in NAB2_fusion):
                NAB2_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in ALK_fusion and i["five_prime_gene"] in ALK_fusion):
                ALK_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])
            
            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in COL1A1_fusion and i["five_prime_gene"] in COL1A1_fusion):
                COL1A1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in ETV6_fusion and i["five_prime_gene"] in ETV6_fusion):
                ETV6_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in FUS_fusion and i["five_prime_gene"] in FUS_fusion):
                FUS_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in EWSR1_fusion and i["five_prime_gene"] in EWSR1_fusion):
                EWSR1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in ATF1_fusion and i["five_prime_gene"] in ATF1_fusion):
                ATF1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])
    
            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in CSF1_fusion or i["five_prime_gene"] in CSF1_fusion):
                CSF1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in WWTR1_fusion and i["five_prime_gene"] in WWTR1_fusion):
                WWTR1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "Cnv" and i["gene_symbol"] in MYC_cnv:
                MYC_result = ["MYC"]

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in PAX3_fusion and i["five_prime_gene"] in PAX3_fusion):
                PAX3_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in NCOA2_fusion and i["five_prime_gene"] in NCOA2_fusion):
                NCOA2_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in HEY1_fusion and i["five_prime_gene"] in HEY1_fusion):
                HEY1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])
            
            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in ZNF444_fusion and i["five_prime_gene"] in ZNF444_fusion):
                ZNF444_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in SS18_fusion and i["five_prime_gene"] in SS18_fusion):
                SS18_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in ASPSCR1_fusion and i["five_prime_gene"] in ASPSCR1_fusion):
                ASPSCR1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])  

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in CREB1_fusion and i["five_prime_gene"] in CREB1_fusion):
                CREB1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in NR4A3_fusion and i["five_prime_gene"] in NR4A3_fusion):
                NR4A3_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"]) 

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in WT1_fusion and i["five_prime_gene"] in WT1_fusion):
                WT1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in BCOR_fusion and i["five_prime_gene"] in BCOR_fusion):
                BCOR_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in TFE3_fusion or i["five_prime_gene"] in TFE3_fusion):
                TFE3_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in ETV1_fusion and i["five_prime_gene"] in ETV1_fusion):
                ETV1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in CIC_fusion and i["five_prime_gene"] in CIC_fusion):
                CIC_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in CCNB3_fusion and i["five_prime_gene"] in CCNB3_fusion):
                CCNB3_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in YWHAE_fusion and i["five_prime_gene"] in YWHAE_fusion):
                YWHAE_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in SP3_fusion and i["five_prime_gene"] in SP3_fusion):
                SP3_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in SQSTM1_fusion and i["five_prime_gene"] in SQSTM1_fusion):
                SQSTM1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in TFCP2_fusion and i["five_prime_gene"] in TFCP2_fusion):
                TFCP2_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in GLI1_fusion and i["five_prime_gene"] in GLI1_fusion):
                GLI1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in MALAT1_fusion and i["five_prime_gene"] in MALAT1_fusion):
                MALAT1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in SRF_fusion and i["five_prime_gene"] in SRF_fusion):
                SRF_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in OGT_fusion and i["five_prime_gene"] in OGT_fusion):
                OGT_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in MEIS1_fusion and i["five_prime_gene"] in MEIS1_fusion):
                MEIS1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in CREB2_fusion and i["five_prime_gene"] in CREB2_fusion):
                CREB2_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in MIR143_fusion and i["five_prime_gene"] in MIR143_fusion):
                MIR143_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in RAF1_fusion and i["five_prime_gene"] in RAF1_fusion):
                RAF1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])   

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in PGR_fusion or i["five_prime_gene"] in PGR_fusion):
                PGR_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in PLAG1_fusion or i["five_prime_gene"] in PLAG1_fusion):
                PLAG1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in USP6_fusion and i["five_prime_gene"] in USP6_fusion):
                USP6_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in NUTM2A_fusion and i["five_prime_gene"] in NUTM2A_fusion):
                NUTM2A_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in JAZF1_fusion and i["five_prime_gene"] in JAZF1_fusion):
                JAZF1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in ESR1_fusion or i["five_prime_gene"] in ESR1_fusion):
                ESR1_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "PSeqRnaSv" and (i["three_prime_gene"] in NCOA3_fusion or i["five_prime_gene"] in NCOA3_fusion):
                NCOA3_result.append(i["five_prime_gene"] + "-" + i["three_prime_gene"])

            if i["bio_category"] == "Snvindel" and i["gene_symbol"] in CTNNB1_gene and not re.search("del|ins|fs", i["hgvs_p"]):
                CTNNB1_result.append("CTNNB1")  

            #if i["bio_category"] == "Snvindel" and i["gene_symbol"] in KRAS_gene and not re.search("del|ins|fs", i["hgvs_p"]):
            #    KRAS_result.append(i["gene_symbol"]) 

            #if i["bio_category"] == "Snvindel" and i["gene_symbol"] in TP53_gene and not re.search("del|ins|fs", i["hgvs_p"]):
            #    TP53_result.append("TP53")

            if i["bio_category"] == "Snvindel" and i["gene_symbol"] in IDH1_gene and not re.search("del|ins|fs", i["hgvs_p"]):
                IDH1_result.append(i["gene_symbol"])

        def dedup(li):
            new = []
            for var in li:
                if var not in new:
                    new.append(var)
            return new   

        self.data_js['result']['NTRK_result'] = dedup(NTRK_result)
        self.data_js['result']['CDK4_result'] = dedup(CDK4_result)
        self.data_js['result']['DDIT3_result'] = dedup(DDIT3_result)
        self.data_js['result']['NAB2_result'] = dedup(NAB2_result)
        self.data_js['result']['ALK_result'] = dedup(ALK_result)
        self.data_js['result']['COL1A1_result'] = dedup(COL1A1_result)
        self.data_js['result']['ETV6_result'] = dedup(ETV6_result)
        self.data_js['result']['FUS_result'] = dedup(FUS_result)
        self.data_js['result']['EWSR1_result'] = dedup(EWSR1_result)
        self.data_js['result']['ATF1_result'] = dedup(ATF1_result)
        self.data_js['result']['CSF1_result'] = dedup(CSF1_result)
        self.data_js['result']['WWTR1_result'] = dedup(WWTR1_result)
        self.data_js['result']['MYC_result'] = dedup(MYC_result)
        self.data_js['result']['PAX3_result'] = dedup(PAX3_result)
        self.data_js['result']['NCOA2_result'] = dedup(NCOA2_result)
        self.data_js['result']['HEY1_result'] = dedup(HEY1_result)
        self.data_js['result']['ZNF444_result'] = dedup(ZNF444_result)
        self.data_js['result']['SS18_result'] = dedup(SS18_result)
        self.data_js['result']['ASPSCR1_result'] = dedup(ASPSCR1_result)
        self.data_js['result']['CREB1_result'] = dedup(CREB1_result)
        self.data_js['result']['WT1_result'] = dedup(WT1_result)
        self.data_js['result']['NR4A3_result'] = dedup(NR4A3_result)
        self.data_js['result']['BCOR_result'] = dedup(BCOR_result)
        self.data_js['result']['TFE3_result'] = dedup(TFE3_result)
        self.data_js['result']['ETV1_result'] = dedup(ETV1_result)
        self.data_js['result']['CIC_result'] = dedup(CIC_result)
        self.data_js['result']['CCNB3_result'] = dedup(CCNB3_result)
        self.data_js['result']['YWHAE_result'] = dedup(YWHAE_result)
        self.data_js['result']['SP3_result'] = dedup(SP3_result)
        self.data_js['result']['SQSTM1_result'] = dedup(SQSTM1_result)
        self.data_js['result']['TFCP2_result'] = dedup(TFCP2_result)
        self.data_js['result']['GLI1_result'] = dedup(GLI1_result)
        self.data_js['result']['MALAT1_result'] = dedup(MALAT1_result)
        self.data_js['result']['SRF_result'] = dedup(SRF_result)
        self.data_js['result']['OGT_result'] = dedup(OGT_result)
        self.data_js['result']['MEIS1_result'] = dedup(MEIS1_result)
        self.data_js['result']['CREB2_result'] = dedup(CREB2_result)
        self.data_js['result']['MIR143_result'] = dedup(MIR143_result)
        self.data_js['result']['RAF1_result'] = dedup(RAF1_result)
        self.data_js['result']['PGR_result'] = dedup(PGR_result)
        self.data_js['result']['PLAG1_result'] = dedup(PLAG1_result)
        self.data_js['result']['USP6_result'] = dedup(USP6_result)
        self.data_js['result']['NUTM2A_result'] = dedup(NUTM2A_result)
        self.data_js['result']['JAZF1_result'] = dedup(JAZF1_result)
        self.data_js['result']['ESR1_result'] = dedup(ESR1_result)
        self.data_js['result']['NCOA3_result'] = dedup(NCOA3_result)
        self.data_js['result']['CTNNB1_result'] = dedup(CTNNB1_result)
        #self.data_js['result']['KRAS_result'] = dedup(KRAS_result)
        #self.data_js['result']['TP53_result'] = dedup(TP53_result)
        self.data_js['result']['IDH1_result'] = dedup(IDH1_result)
        return 

    def get_msi_info(self):
        """
        msi
        :return:
        """
        # msi
        data = self.data_js
        msi = data['msi']
        if not msi:  # 如果为空
            self.data_js['result']['msi'] = []
            return
        for i in msi:
            i['gene_symbol'] = 'MSI'
            i['hgvs_p'] = i['var_id']
            i['class'] = '-'
        if msi[0]['hgvs_p'] == 'MSS':
            msi = [{'var_id' : 'MSS'}]
        else:
            msi[0]['regimen_drug'] = [i for i in msi[0]['evi_sum'] if i['regimen_name']]

        self.data_js['result']['msi'] = msi[0] if msi else {}

        return self.data_js['result']['msi']

    def get_qc_info(self):
        """
        QC
        :return:
        """
        # ===== DNA 质控 =====
        # 湿实验质控
        lib_dna_qc = self.data_js['lib_quality_control']['lib_dna_qc'] if type(self.data_js['lib_quality_control']['lib_dna_qc']).__name__ == 'dict' else self.data_js['lib_quality_control']['lib_dna_qc'][0]
        # 简化，不用self.data_js['lib_quality_control']['lib_dna_qc']判断
        #total_dna = lib_dna_qc.get('dna_qty', '') if self.data_js['lib_quality_control']['lib_dna_qc'] else ''
        #total_dna_lib = lib_dna_qc.get('library_qty', '') if self.data_js['lib_quality_control']['lib_dna_qc'] else ''
        #dna_frag_len = lib_dna_qc.get('library_fragment_length', '') if self.data_js['lib_quality_control']['lib_dna_qc'] else ''

        total_dna = lib_dna_qc.get('dna_qty', '') if lib_dna_qc else ''
        total_dna_lib = lib_dna_qc.get('library_qty', '') if lib_dna_qc else ''
        dna_frag_len = lib_dna_qc.get('library_fragment_length', '') if lib_dna_qc else ''

        # 干实验质控
        dna_qc = self.data_js['qc']['dna_data_qc'] if type(self.data_js['qc']['dna_data_qc']).__name__ == 'dict' \
            else self.data_js['qc']['dna_data_qc'][0] if self.data_js['qc']['dna_data_qc'] else []
        dna_depth_ssbc = dna_qc.get('depth_ssbc', '') if dna_qc else ''
        dna_uni20 = dna_qc.get('uni20', '') if dna_qc else ''
        dna_mapping_ratio = dna_qc.get('mapping_ratio', '') if dna_qc else ''
        dna_cover_ratio = dna_qc.get('cover_ratio', '') if dna_qc else ''
        
        dna_q30 = dna_qc.get('cleandata_q30', '') if dna_qc else ''

        # ===== RNA 质控 =====
        # 湿实验质控
        lib_rna_qc = self.data_js['lib_quality_control']['rna_lib_qc'] if type(self.data_js['lib_quality_control']['rna_lib_qc']).__name__ == 'dict' else self.data_js['lib_quality_control']['rna_lib_qc'][0]
        # 两行填表，应该从rna_lib_qc抓取质控数值
        #total_rna = lib_dna_qc.get('rna_qty', '') if self.data_js['lib_quality_control']['lib_dna_qc'] else ''
        # 填RNA文库总量
        #total_rna_lib = lib_dna_qc.get('rna_fnl_library_qty', '') if self.data_js['lib_quality_control']['lib_dna_qc'] else ''
        # 要求填RNA预文库总量，2023年7月14日
        #total_rna_lib = lib_dna_qc.get('rna_pre_library_qty', '') if self.data_js['lib_quality_control']['lib_dna_qc'] else ''

        total_rna = lib_rna_qc.get('rna_qty', '') if lib_rna_qc else ''
        total_rna_lib = lib_rna_qc.get('rna_pre_library_qty', '') if lib_rna_qc else ''
        # 不做RNA文库片段长度质控，也用不到，不删了
        #rna_frag_len = lib_dna_qc.get('rna_fnl_library_fragment_length', '') if self.data_js['lib_quality_control']['lib_dna_qc'] else ''
        
        # 干实验质控
        rna_qc = self.data_js['qc']['rna_data_qc'] if type(self.data_js['qc']['rna_data_qc']).__name__ == 'dict' \
            else self.data_js['qc']['rna_data_qc'][0] if self.data_js['qc']['rna_data_qc'] else []
        rna_q30 = rna_qc.get('cleandata_q30', '') if rna_qc else ''
        rna_cleandata = rna_qc.get('cleandata_size', '') if rna_qc else ''
        rrna_ratio = rna_qc.get('rrna_ratio', '') if rna_qc else ''

        # 方便核对报告模板中的质控抓取字段，重命名湿实验相关key
        self.data_js['result']['qc'] = {
            'dna_qty' : '%.2f' % (round(float(total_dna), 2)) if total_dna else '',
            'library_qty' : '%.2f' % (round(float(total_dna_lib), 2)) if total_dna_lib else '',
            'library_fragment_length' : dna_frag_len if dna_frag_len else '',
            'dna_cover_ratio': '%.2f' % (round(float(dna_cover_ratio) * 100, 2)) if dna_cover_ratio else '',
            #'dna_q30': '%.2f' % (round(float(dna_q30 * 100), 2)) if dna_q30 else '',
            'dna_q30' : '%.2f' % (round(float(dna_q30) * 100, 2)) if dna_q30 else '',
            'dna_depth_ssbc': '%.2f' % round(float(dna_depth_ssbc), 2) if dna_depth_ssbc else '',
            'dna_uni20': '%.2f' % round(float(dna_uni20) * 100, 2) if dna_uni20 else '',
            'dna_mapping_ratio': '%.2f' % round(float(dna_mapping_ratio) * 100, 2) if dna_mapping_ratio else '',
            # RNA
            'rna_q30': '%.2f' % (round(float(rna_q30) * 100, 2)) if rna_q30 else '',          
            'rna_qty': '%.2f' % (round(float(total_rna), 2)) if total_rna else '',
            'rna_pre_library_qty': '%.2f' % (round(float(total_rna_lib), 2)) if total_rna_lib else '',
            'rna_cleandata': '%.2f' % float(rna_cleandata[:-1]) if rna_cleandata else '',
            'rrna_ratio': '%.2f' % round(float(rrna_ratio) * 100, 2) if rrna_ratio else '',
        }

    def get_reference(self):
        """
        参考文献
        :return:
        """
        self.data_js['result']['reference'] = self.data_js['refer']

        return self.data_js['result']['reference']

    def to_word(self):
        """
        写入模板
        :return:
        """
        base_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(base_dir, 'template', 'RNASeqSarcoma-CustomEdition-hospital-ZJZL-v14.docx')
        tpl = docxtpl.DocxTemplate(path)
        tpl.render(self.data_js)
        tpl.save(os.path.join(self.output_div, self.json_name + '.docx'))
        tpl.save(os.path.join(self.output_div, self.json_name + '_RZL.docx'))
        
        return 
    
    
    def to_json(self):
        '''
        在json文件中追加变异统计信息
        '''
        msi = self.get_msi_info()
        if msi:
            msi_content = '微卫星稳定型（MSS）' if msi['var_id'] == 'MSS' else '微卫星不稳定型（MSI-H）'
        else:
            msi_content = ''

        self.data_js['var_summary'] = {}
        var_sum = self._var_count()

        fgene_content = '检出' + str(var_sum['somatic_all'] + var_sum['germline_4_5']) + '个基因变异：' + '\n' +'其中靶向用药相关的变异有' + str(var_sum['predictive']) + '个；' + '\n'\
            + '诊断相关的变异有' + str(var_sum['diagnostic']) + '个' if var_sum['somatic_all'] + var_sum['germline_4_5'] else '未检出相关变异'
        fsomatic_cell = '检出' + str(var_sum['somatic_all']) + '个基因变异：' + '\n' +'其中靶向用药相关的变异有' + str(var_sum['predictive']) + '个;' + '\n' + '诊断相关的变异有' + str(var_sum['diagnostic']) + '个' if var_sum['somatic_all'] else '未检出相关变异'
        fblast_cell = '检出' + str(var_sum['germline_4_5']) + '个致病变异或疑似致病变异，其中' + str(var_sum['germline_drug']) + '个与靶向药物相关' if var_sum['germline_4_5'] else '未检出相关变异'
        var_summary = {'fgene_content' : fgene_content,
                        'fsomatic_cell' : fsomatic_cell,
                        'fblast_cell' : fblast_cell,
                        'fmsi_content' : msi_content,
                        'ftmb_content' : '',
                        'report_remark' : ''}

        self.data_js['var_summary'] = var_summary
        dataJson = json.dumps(self.data_js, ensure_ascii=False)
        with open(self.output_div + '/' + self.json_name + '_for_Sar_summary.json', 'w', encoding='utf-8') as outFile:
            outFile.write(dataJson)

    

    def run(self):
        """
        程序入口
        :return:
        """
        # 基本信息
        self.get_base_info()
        # # 检测结果 (检测小结)
        self.get_test_result()
        # # 基因变异结果汇总
        self.get_var()
        # MSI
        self.get_msi_info()
        # # qc
        self.get_qc_info()
        # # 参考文献
        self.get_reference()
        # 肉瘤
        #self.get_rouliu_result()
        self.get_sarcoma_result()
        # 写入模板
        self.to_word()
        self.to_json()

        return


def parse_args():
    """
    outfile是一个目录, 存放json的目录
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--json_name', dest='json_name', required=True)
    parser.add_argument('-o', '--outfile', dest='outfile', required=True)
    arg = parser.parse_args()
    return arg


if __name__ == '__main__':     
    args = parse_args()
    report = Base(args.json_name, output_div=args.outfile)
    report.run()