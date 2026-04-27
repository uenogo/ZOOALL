#coding: UTF-8
"""
(C) RIKEN/JASRI 2020 :
Author: Kunio Hirata
ESA -> the main function : the class to read & write zoo database file
This code is originally written by K.Hirata and modified by N.Mizuno.
NM added function to read xlsx file directly and output zoo.db by using ESA class.

Vice-author: Nobuhiro Mizuno
"""
import sys, os, math, numpy, csv, re, datetime, xlrd, codecs
import configparser
import pandas as pd
import numpy as np
import KUMA
# logger の設定
import logging
from configparser import ConfigParser, ExtendedInterpolation

class UserESA():
    def __init__(self, fname=None, root_dir=".", beamline=None):
        # beamlineの名前はconfigから読む
        self.config = ConfigParser(interpolation=ExtendedInterpolation())
        config_path = "%s/beamline.ini" % os.environ['ZOOCONFIGPATH']
        self.config.read(config_path)

        self.fname = fname
        self.isRead = None
        self.isPrep = None
        self.isGot  = None
        self.zoocsv = None
        self.contents = []

        # configure file から情報を読む: beamlineの名前
        self.beamline = self.config.get("beamline", "beamline")
        import BeamsizeConfig
        self.bsconf = BeamsizeConfig.BeamsizeConfig()

        # CSV prefix
        self.csv_prefix = self.fname.replace(".xlsx","")

        # logger の設定
        self.logger = logging.getLogger("ZOO")
        # output log string with level 'info'
        self.logger.setLevel(logging.INFO)
        # levelがwarningのときには標準出力とファイル両方に出力する
        
        # create file handler which logs even debug messages
        self.logger_fh = logging.FileHandler('useresa.log')
        self.logger_fh.setLevel(logging.DEBUG)
        # create console handler with a higher log level
        self.logger_ch = logging.StreamHandler()
        self.logger_ch.setLevel(logging.WARNING)
        # create formatter and add it to the handlers
        self.logger_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # set formatter to handlers
        self.logger_fh.setFormatter(self.logger_formatter)
        self.logger_ch.setFormatter(self.logger_formatter)
        # add the handlers to logger
        self.logger.addHandler(self.logger_fh)
        self.logger.addHandler(self.logger_ch)

        self.root_dir = root_dir

    def setDefaults(self):
        # self.df に以下のカラムを追加する
        # self.config.getfloat("experiment", "score_min") などで読み込む
        # "score_min"
        # "score_max"
        # "raster_dose"
        # "dose_ds"
        # "raster_roi"
        # "exp_raster"
        # "att_raster"
        # "hebi_att"
        # "cover_flag"
        # "exp_ds"
        self.df["score_min"] = self.config.getfloat("experiment", "score_min")
        self.df["score_max"] = self.config.getfloat("experiment", "score_max")
        self.df["raster_dose"] = self.config.getfloat("experiment", "raster_dose")
        self.df["dose_ds"] = self.config.getfloat("experiment", "dose_ds")
        self.df["raster_roi"] = self.config.getint("experiment", "raster_roi")
        self.df["exp_ds"] = self.config.getfloat("experiment", "exp_ds")
        self.df["exp_raster"] = self.config.getfloat("experiment", "exp_raster")
        # att_raster の数値を取得して小数点以下第一位までに丸める
        self.df["att_raster"] = self.config.getfloat("experiment", "att_raster")
        self.df["att_raster"] = round(self.df["att_raster"], 1)
        self.df["hebi_att"] = self.df["att_raster"]
        self.df["cover_scan_flag"] = self.config.getint("experiment", "cover_flag")
        # 結晶サイズは max_crystal_size として読み込む
        self.df['cry_min_size_um'] = self.df['max_crystal_size']
        self.df['cry_max_size_um'] = self.df['max_crystal_size']
        # root_dir は self.root_dir として読み込む
        self.df['root_dir'] = self.root_dir
        # p_indexはDataFrameのインデックスと同じで良い
        self.df['p_index'] = self.df.index
        # offset_angle は 0 とする
        self.df['offset_angle'] = 0
        # reduced_fact は 1 とする
        self.df['reduced_fact'] = 1
        # ntimes は 1 とする
        self.df['ntimes'] = 1
        # meas_name は 変換に利用したファイル名を入れておく
        self.df['meas_name'] = self.fname 
        # hel_full_osc,hel_part_osc
        self.df['hel_full_osc'] = 60.0
        self.df['hel_part_osc'] = 30.0

        # 'desired_exp' と 'mode' から実験パラメータを設定する
        # 1) desired_exp が "scan_only" のとき 
        # score_min, score_max ともに 9999 とする
        # raster_dose: 0.3, dose_ds: 0.0, cover_flag: 0
        self.df.loc[self.df['desired_exp'] == "scan_only", 'score_min'] = 9999
        self.df.loc[self.df['desired_exp'] == "scan_only", 'raster_dose'] = 0.3
        self.df.loc[self.df['desired_exp'] == "scan_only", 'dose_ds'] = 0.0
        self.df.loc[self.df['desired_exp'] == "scan_only", 'cover_scan_flag'] = 0

        # 2) desired_exp が "normal" のとき
        # mode が "helical" の場合には、score_max を 9999 とする
        self.df.loc[self.df['desired_exp'] == "normal", 'score_max'] = 9999

        # 3) desired_exp が "ultra_high_dose_scan" のとき
        # dose_ds = 9.0 とする
        self.df.loc[self.df['desired_exp'] == "ultra_high_dose_scan", 'dose_ds'] = 9.0
        # 4) desired_exp が "phaing" のとき
        # dose_ds = 5.0 とする
        self.df.loc[self.df['desired_exp'] == "phasing", 'dose_ds'] = 5.0

    # ビームライン、実験モードと結晶のタイプから実験パラメータを取得する
    # 2023/05/09 type_crystal は使わない
    def getParams(self, desired_exp_string, mode):
        desired_exp_string = desired_exp_string.lower()

        # DEFAULT PARAMETER
        # beamline.ini から読む
        #self.beamline = self.config.get("beamline", "beamline")
        score_min   = self.config.getfloat("experiment", "score_min")
        score_max   = self.config.getfloat("experiment", "score_max")
        raster_dose = self.config.getfloat("experiment", "raster_dose")
        dose_ds     = self.config.getfloat("experiment", "dose_ds")
        raster_roi  = self.config.getint("experiment", "raster_roi")
        exp_raster = self.config.getfloat("experiment", "exp_raster")
        att_raster  = self.config.getfloat("experiment", "att_raster")
        hebi_att    = self.config.getfloat("experiment", "hebi_att")
        cover_flag  = self.config.getint("experiment", "cover_flag")

        # PARAMTER CONDITION
        self.param = {
            "scan_only":{
                "single":   [9999, 9999, 0.3, dose_ds, 0, exp_raster, att_raster, hebi_att, 0],
                "helical":  [9999, 9999, 0.3, dose_ds, 0, exp_raster, att_raster, hebi_att, 0],
                "multi":    [9999, 9999, 0.3, dose_ds, 0, exp_raster, att_raster, hebi_att, 0],
                "mixed":    [9999, 9999, 0.3, dose_ds, 0, exp_raster, att_raster, hebi_att, 0],
            },
            "normal":{
                "single":   [score_min, score_max, 0.1, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "helical":  [score_min, 9999, 0.05, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "multi":    [score_min, score_max, 0.1, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "mixed":    [score_min, 9999, 0.1, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
            },
            "high_dose_scan":{
                "single":   [score_min, 9999, 0.05, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "helical":  [score_min, 9999, 0.05, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "multi":    [score_min, 9999, 0.05, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "mixed":    [score_min, 9999, 0.05, dose_ds, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
            },
            "ultra_high_dose_scan":{
                "single":   [score_min, score_max, 0.2, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
                "helical":  [score_min, score_max, 0.2, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
                "multi":    [score_min, score_max, 0.2, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
                "mixed":    [score_min, score_max, 0.2, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
            },
            "phasing":{
                "single":   [score_min, score_max, 0.1, 5, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "helical":  [score_min, 9999, 0.05, 5, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "multi":    [score_min, score_max, 0.1, 5, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
                "mixed":    [score_min, score_max, 0.1, 5, raster_roi, exp_raster, att_raster, hebi_att, cover_flag],
            },
            "rapid":{
                "single":   [score_min, score_max, raster_dose, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
                "helical":  [score_min, score_max, raster_dose, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
                "multi":    [score_min, score_max, raster_dose, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
                "mixed":    [score_min, score_max, raster_dose, dose_ds, raster_roi, exp_raster, 100, 100, cover_flag],
            },
        }

        return self.param[desired_exp_string][mode]
    
    def checkLN2flag(self):
        # self.dfのカラム "ln2_flag" について以下のパターンで処理を行う
        # 'NaN'であれば ０
        # 'Yes' or 'yes' or "YES" であれば １
        # 'Unavailable' であれば ０
        # それ以外であれば ０
        self.df['ln2_flag'] = self.df['ln2_flag'].fillna(0)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('Yes', 1)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('yes', 1)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('YES', 1)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('Unavailable', 0)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('No', 0)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('no', 0)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('NO', 0)
        self.df['ln2_flag'] = self.df['ln2_flag'].replace('-', 0)

        #print(self.df)

    def checkZoomFlag(self):
        # self.dfのカラム "ln2_flag" について以下のパターンで処理を行う
        # 'NaN'であれば ０
        # 'Yes' or 'yes' or "YES" であれば １
        # 'Unavailable' であれば ０
        # それ以外であれば ０
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].fillna(0)
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].replace('Yes', 1)
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].replace('yes', 1)
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].replace('YES', 1)
        # self.df['zoomcap_flag']が　'No' or 'no' or 'NO' or 'Unavailable' であれば ０
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].replace('No', 0)
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].replace('no', 0)
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].replace('NO', 0)
        self.df['zoomcap_flag'] = self.df['zoomcap_flag'].replace('Unavailable', 0)

        # DataFrameを省略することなく表示する
        pd.set_option('display.max_rows', None)
        #print(self.df)

    def checkPinFlag(self):
        #print(self.df['pin_flag'])
        # self.df['warm_time']の初期値を30.0とする
        self.df['warm_time'] = 30.0
        # self.df にはすでに"pin_flag"があるので、それを利用する
        # self.df['pin_flag']の文字列を小文字に変換した文字列が "spine"　であれば self.df['warm_time'] = 10.0
        self.df.loc[self.df['pin_flag'].str.lower() == 'spine', 'warm_time'] = 10.0
        # self.df['pin_flag']の文字列を小文字に変換した文字列が "als + ssrl"　であれば self.df['warm_time'] = 20.0
        self.df.loc[self.df['pin_flag'].str.lower() == 'als + ssrl', 'warm_time'] = 20.0
        # self.df['pin_flag']の文字列を小文字に変換した文字列が "copper"　であれば self.df['warm_time'] = 60.0
        self.df.loc[self.df['pin_flag'].str.lower() == 'copper', 'warm_time'] = 60.0
        # self.df['pin_flag']の文字列を小文字に変換した文字列が "no-wait"　であれば self.df['warm_time'] = 0.0
        self.df.loc[self.df['pin_flag'].str.lower() == 'no-wait', 'warm_time'] = 0.0

    def fillFlux(self):
        # self.df['flux']の数値を読み込む
        # self.bsconf.getFluxAtWavelength(hbeam, vbeam, wavelength)を呼び出す
        # この関数の引数に self.df['hbeam'], self.df['vbeam'], self.df['wavelength']を渡す
        # 戻り値はfluxである
        # fluxの値をself.df['flux']に代入する
        self.df['flux'] = self.df.apply(lambda x: self.bsconf.getFluxAtWavelength(x['ds_hbeam'], x['ds_vbeam'], x['wavelength']), axis=1)

    def splitBeamsizeInfo(self):
        # self.df['beamsize']の文字列をself.checkBeamsizeの引数として渡す
        # self.checkBeamsize()は self.df['beamsize']を引数とし、戻り値は(hbeam, vbeam)である(どちらもfloatのタプル)
        # hbeam, vbeamの数値は新たなカラムとしてself.dfに追加される 'hbeam', 'vbeam'
        self.df['ds_hbeam'], self.df['ds_vbeam'] = zip(*self.df['beamsize'].map(self.checkBeamsize))
        self.df['raster_hbeam'], self.df['raster_vbeam'] = zip(*self.df['beamsize'].map(self.checkBeamsize))

    # Raster scanの露光条件を定義する
    # Pandas dataframeに対して一気に処理を行う
    def defineScanCondition(self):
        # Dose estimation will be conducted by KUMA
        kuma = KUMA.KUMA()
    
        # self.df['wavelgnth']からself.df['energy']を計算する
        self.df['energy'] = 12.3984 / self.df['wavelength']

        # self.df['desired_exp']の文字列を小文字に変換した文字列が "normal", "scan_only", "phasing", "rapid"の場合は以下の処理を行う
        # photons_per_image は 4E10 で固定する
        # photons_per_exptime は flux * exp_raster で計算する
        # df['att_raster'] = photons_per_image / photons_per_exptime * 100.0 とする
        # df['hebi_att'] = photons_per_image / photons_per_exptime * 100.0 とする
        # maskを利用して条件ごとに処理をしていく
        mask1 = (self.df['desired_exp'] == 'normal') | (self.df['desired_exp'] == 'scan_only') | (self.df['desired_exp'] == 'phasing') | (self.df['desired_exp'] == 'rapid')
        photons_per_image = 4E10
        # 1 secあたりの最大フォトン数を計算する
        photons_per_exptime = self.df['flux'] * self.df['exp_raster']
        # 1 frameあたりに必要なフォトン数を入れるためのatt_factorを計算する
        self.df.loc[mask1, 'att_raster'] = photons_per_image / photons_per_exptime * 100.0
        # 1 frameあたりに必要なフォトン数を入れるためのhebi_att_factorを計算する
        self.df.loc[mask1, 'hebi_att'] = photons_per_image / photons_per_exptime * 100.0
        # 1 frameあたりのphotonsを計算する
        self.df.loc[mask1, 'ppf_raster'] = photons_per_image
        # 1 frameあたりのdoseを計算する
        # kuma.getDose()の引数は hbeam, vbeam, flux, energy, exp_raster
        # dose_per_frame = kuma.getDose(hbeam, vbeam, flux, energy, exp_raster) * self.df['att_raster'] / 100.0
        self.df.loc[mask1, 'dose_per_frame'] = kuma.getDose(self.df['ds_hbeam'], self.df['ds_vbeam'], self.df['flux'], self.df['energy'], self.df['exp_raster']) * self.df['att_raster'] / 100.0

        # mask2 
        mask2 = (self.df['desired_exp'] == 'high_dose_scan')
        dose_for_raster = 0.30 # MGy
        # 1 frameあたりのdoseを計算する
        self.df.loc[mask2, 'dose_per_frame'] = kuma.getDose(self.df['ds_hbeam'], self.df['ds_vbeam'], self.df['flux'], self.df['energy'], self.df['exp_raster'])
        # transmissionは dose_for_raster / dose_per_frame * 100.0 で計算する
        self.df.loc[mask2, 'att_raster'] = dose_for_raster / self.df['dose_per_frame'] * 100.0
        self.df.loc[mask2, 'hebi_att'] = dose_for_raster / self.df['dose_per_frame'] * 100.0
        # 'ppf' = photons per frame
        self.df.loc[mask2, 'ppf_raster'] = self.df['flux'] * self.df['exp_raster'] * self.df['att_raster'] / 100.0
        # dose_per_frame = kuma.getDose(hbeam, vbeam, flux, energy, exp_raster) * self.df['att_raster'] / 100.0
        self.df.loc[mask2, 'dose_per_frame'] = kuma.getDose(self.df['ds_hbeam'], self.df['ds_vbeam'], self.df['flux'], self.df['energy'], self.df['exp_raster']) * self.df['att_raster'] / 100.0

        # masks
        mask3 = (self.df['desired_exp'] == 'ultra_high_dose_scan')
        dose_for_raster = 1.0 # MGy
        # 1 frame あたりのdoseを計算する
        self.df.loc[mask3, 'dose_per_frame'] = kuma.getDose(self.df['ds_hbeam'], self.df['ds_vbeam'], self.df['flux'], self.df['energy'], self.df['exp_raster'])
        # transmissionは dose_for_raster / dose_per_frame * 100.0 で計算する
        self.df.loc[mask3, 'att_raster'] = dose_for_raster / self.df['dose_per_frame'] * 100.0
        self.df.loc[mask3, 'hebi_att'] = dose_for_raster / self.df['dose_per_frame'] * 100.0
        # 'ppf' = photons per frame
        self.df.loc[mask3, 'ppf_raster'] = self.df['flux'] * self.df['exp_raster'] * self.df['att_raster'] / 100.0
        # dose_per_frame = kuma.getDose(hbeam, vbeam, flux, energy, exp_raster) * self.df['att_raster'] / 100.0
        self.df.loc[mask3, 'dose_per_frame'] = kuma.getDose(self.df['ds_hbeam'], self.df['ds_vbeam'], self.df['flux'], self.df['energy'], self.df['exp_raster']) * self.df['att_raster'] / 100.0

        # self.logger.info -> 'dose_per_frame' をリスト表示
        # puckid, pinid, dose_per_frame のリストを表示する
        # format f"PuckID: {puckid} PinID: {pinid} dose_per_frame: {dose_per_frame}"
        # 'pinid' については文字列の場合があるのでそのまま表示する
        self.logger.info("Scan conditions estimated results")
        for i in range(len(self.df)):
            self.logger.info(f"PuckID: {self.df['puckid'][i]} PinID: {self.df['pinid'][i]} dose_per_frame: {self.df['dose_per_frame'][i]:.3f}")
            #self.logger.info(f"Puck-Pin ID: {self.df['puckid'][i]}-{self.df['pinid'][i]:2d} : dose/frame: {self.df['dose_per_frame'][i]:.3f} MGy")

        #print(self.df)


    # end of defineScanCondition()

    def makeExpWarning(self): 
        # 1 frameあたりのdoseが0.3MGyを超えていて、self.df['desired_exp'] が 'high_dose_scan' もしくは 'ultra_high_dose_scan'出ない場合は警告を出す
        # 丁寧な文字列でloggerを出力する
        # "Warning: dose/frame exceeds 0.3 MGy. Please check the exposure condition."
        # "puckid: 'sample' pinid: 01 dose/frame 0.5 MGy"
        mask = (self.df['dose_per_frame'] > 0.3) & (self.df['desired_exp'] != 'high_dose_scan') & (self.df['desired_exp'] != 'ultra_high_dose_scan')
        
        if mask.any():
            for i in range(len(self.df)):
                if mask[i]:
                    self.logger.warning("Warning: dose/frame exceeds 0.3 MGy. Please check the exposure condition.")
                    self.logger.warning("puckid: {} pinid: {} dose/frame {} MGy".format(self.df['puckid'][i], self.df['pinid'][i], self.df['dose_per_frame'][i]))
        else:
            self.logger.info("No warning message for dose/frame check.")

        # self.df['ppf_raster']が 4.0E10 を下回る場合には警告を出す
        # "Warning: ppf_raster is less than 4.0E10. Please check the exposure condition."
        mask2 = (self.df['ppf_raster'] < 4.0E10)
        if mask2.any():
            for i in range(len(self.df)):
                if mask2[i]:
                    self.logger.warning("Warning: ppf_raster is less than 4.0E10. Please check the exposure condition.")
                    self.logger.warning("puckid: {} pinid: {} ppf_raster {}".format(self.df['puckid'][i], self.df['pinid'][i], self.df['ppf_raster'][i]))
        else:
            self.logger.info("No warning message for photons/frame check")

    def sizeWarning(self):
        # self.df['mode']が 'multi' である場合、self.df['max_crystal_size']と self.df['beam_size']を比較して、
        # self.df['hbeam']と self.df['vbeam']を比較して大きい方を tmp_beamsize とする
        # self.df['max_crystal_size']が tmp_beamsize の2倍よりも大きい場合には警告を出す
        # "Warning: max_crystal_size is larger than 2 times of beam_size. Please check the exposure condition."
        mask = (self.df['mode'] == 'multi')
        if mask.any():
            for i in range(len(self.df)):
                if mask[i]:
                    tmp_beamsize = max(self.df['ds_hbeam'][i], self.df['ds_vbeam'][i])
                    if self.df['max_crystal_size'][i] > tmp_beamsize * 2.0:
                        self.logger.warning("Warning: max_crystal_size is larger than 2 times of the larger dimension of the beam size.")
                        self.logger.warning("Please re-confirm the conditions of 'multi' mode.")
                        self.logger.warning("puckid: {} pinid: {} max_crystal_size:{:.1f}um beam size:{}um".format(self.df['puckid'][i], self.df['pinid'][i], self.df['max_crystal_size'][i], tmp_beamsize))
        
    # self.dfに格納されているから、データexp_rasterに変更を加える必要がある場合には変更を加える
    def modifyExposureConditions(self):
        # self.df['att_raster']　が 100.0 を超えている場合
        # さらにself.df['exp_raster']を長くして、その分 self.df['att_raster'] = 100.0とする
        # その場合、self.df['hebi_att']も変更する必要がある
        # extend_ratio = self.df['att_raster'] / 100.0
        # new_exp_raster = self.df['exp_raster'] * extend_ratio
        # この数値を self.df['exp_raster'] に代入する
        mask = (self.df['att_raster'] > 100.0)
        self.df.loc[mask, 'exp_raster'] = self.df['exp_raster'] * self.df['att_raster'] / 100.0
        self.df.loc[mask, 'att_raster'] = 100.0
        self.df.loc[mask, 'hebi_att'] = 100.0
        # self.loggerにWarningを出す
        # mask が Trueの場合のみ、Warningを出す
        # そのとき 'puckid', 'pinid' を出力する
        # さらに exp_raster の数値も同時に出力する
        if mask.any():
            self.logger.warning("att_raster > 100.0 -> 'exp_raster' was modified")
            self.logger.warning("Please carefully check 'beam size' and 'desired experimental mode'")
            self.logger.warning(self.df.loc[mask, ['puckid', 'pinid', 'sample_name', 'exp_raster']])

        # self.dfに含まれる露光条件で
        # ppf_rasterが 4.0E10 を下回る場合
        # dose_per_frameが 0.3 MGy を超える場合 にWarning messageを出す
        # loggingに記録する
        self.makeExpWarning()

    def makeCSV(self, zoo_csv=None):
        if not zoo_csv:
            return None

        ctime=datetime.datetime.now()
        time_str = datetime.datetime.strftime(ctime, '%y%m%d%H%M%S')
        db_fname = "zoo_%s.db"%time_str

        import ESA
        esa = ESA.ESA(db_fname)
        if os.path.exists(self.csvout):
            esa.makeTable(self.csvout, force_to_make=True)

        return

    def read_new(self):
        # pandasを利用して.xlsxファイルを読み込む
        # tabの名前を指定して読む "ZOOPREP_YYMMDD_NAME_BLNAME_v2"
        # pandasを利用してエクセルのタブのリストを取得して表示する
        #print(pd.ExcelFile(self.fname).sheet_names)

        # エクセルのタブ名が "ZOOPREP_YYMMDD_NAME_BLNAME_v2" であるタブを読み込む
        # Index(['PuckID', 'PinID', 'SampleName', 'Objective', 'Mode', 'HA',
        # 'Wavelength [Å]', 'Hor. scan length [µm]', 'Resolution limit [Å]',
        # 'Beam size [um]\n(H x V)', 'Crystal size [µm]',
        # '# of crystals\n / Loop', 'Total osc \n/ Crystal', 'Osc. Width',
        # 'LN2\nSplash', 'PIN Type', 'Zoom\nCapture', 'Unnamed: 17',
        # 'Confirmation required'],
        # column名を指定する
        columns = ['puckid', 'pinid', 'sample_name', 'desired_exp', 'mode', 'anomalous_flag', \
            'wavelength', 'loopsize', 'resolution_limit', 'beamsize', 'max_crystal_size', 'maxhits', 'total_osc', 'osc_width', \
                'ln2_flag', 'pin_flag', 'zoomcap_flag', 'what', 'confirmation_require']

        # データは4行目から
        # 250121: sheet_name -> ZOO_YYMMDD_NAME_BLNAME_v2 -> Sheet
        self.df = pd.read_excel(self.fname, sheet_name="Sheet", header=2)
        # 列名を指定する
        self.df.columns = columns
        # 現時点でのデータ数をself.loggerに出力する
        self.logger.info("Number of data: %d"%len(self.df))
        # 'puckid' がないデータを削除する
        self.df = self.df.dropna(subset=['puckid'])
        # 現時点でのデータ数をself.loggerに出力する
        self.logger.info("Number of data after polishment: %d"%len(self.df))
        self.isPrep = True

    def expandPinRange(self, pinstr):
        # pinid_str = "1-4" のような文字列を受け取る
        # 1-4 の場合は、1,2,3,4 のリストを返す
        # 1+2+3 の場合は、1,2,3 のリストを返す
        # 1;2;3 の場合は、1,2,3 のリストを返す
        # それ以外はそのままの文字列を返す
        if '-' in pinstr:
            start, end = map(int, pinstr.split('-'))
            return list(range(start, end + 1))
        elif '+' in pinstr:
            return list(map(int, pinstr.split('+')))
        elif ';' in pinstr:
            return list(map(int, pinstr.split(';')))
        else:
            return [pinstr]

    def dividePinInfo(self, pin_char):
        import re
        # ステップ1: 区切り文字で分割
        parts = re.split(r'[;+.]', pin_char)
        
        # ステップ2: 各部分を range に変換
        ranges = []
        for part in parts:
            if '-' in part:
                start, end = map(int, part.split('-'))
                ranges.append(range(start, end + 1))
            else:
                num = int(part)
                ranges.append(range(num, num + 1))
        
        # ステップ3: 必要ならすべての値をフラットなリストに
        flattened = [i for r in ranges for i in r]
        print(f"flattened={flattened}")  # [1, 2, 3, 4, 5, 10, 11, 12, 15, 16]
        return flattened

    def expandCompressedPinInfo(self):
        # The new dataframe of expanded pins
        new_df_list = []
        isFound=False
        for i, row in self.df.iterrows():
            pinid_str = str(row['pinid'])
            print(f"##### pinid_str= {pinid_str}")
            flattened_ids = self.dividePinInfo(pinid_str)
            for pinid in flattened_ids:
                new_row = row.copy()
                new_row['pinid'] = int(pinid)
                new_df_list.append(new_row)

        # 4. Create a new dataframe from the list of expanded rows
        new_df = pd.DataFrame(new_df_list)
        # 5. Reset the index of the new dataframe
        new_df.reset_index(drop=True, inplace=True)
        # 6. new_df -> self.df
        self.df = new_df

    def calcDist(self, wavelength, resolution_limit, isROI=False):
        # beamline.ini　の experiment セクション　から min_camera_lim を読んで min_camera_len に代入する
        min_camera_len = self.config.getfloat("detector", "min_camera_len")

        # ROIがない場合
        if isROI == False:
            self.logger.info(f"ROI is False")
            # wavelength と resolution_limit から camera_len を計算する
            # camera_len が min_camera_len 以下なら min_camera_len を返す
            min_camera_dim = self.config.getfloat("detector", "min_camera_dim")
        else:
            self.logger.info(f"ROI is True")
            # ROIがある場合なんだが、calcDistFromLength()は半径でなく直径を要求する -> min_camera_dim = 2 * min_camera_dim
            min_camera_dim = self.config.getfloat("experiment", "raster_roi_edge_mm") * 2.0

        camera_len = self.calcDistFromLength(wavelength, resolution_limit, min_camera_dim)

        self.logger.info(f"calcuated camera_len: {camera_len}")

        # camera_len が　min_camera_len 以下なら min_camera_len を返す
        # camera_len が min_dim より大きいなら camera_len を返す
        if camera_len < min_camera_len:
            camera_len = min_camera_len

        # 小数点第一位に丸める camera_len
        camera_len = round(camera_len, 1)

        return camera_len

    def calcDistFromLength(self, wavelength, resolution_limit, detector_radius):
        # wavelength と resolution_limit から camera_len を計算する
        # camera_len が min_dim 以下なら min_dim を返す
        # camera_len が min_dim より大きいなら camera_len を返す
        theta = numpy.arcsin(wavelength / 2.0 / resolution_limit)
        bunbo = 2.0 * numpy.tan(2.0 * theta)
        camera_len = detector_radius / bunbo
        return camera_len

    def checkBeamsize(self, beamsize_char):
        cols = beamsize_char.split('x')
        if len(cols) > 1:
            hbeam = float(cols[0])
            vbeam = float(cols[1])
            return hbeam, vbeam

    # データフレームの分解能限界からカメラ長を計算して格納する
    def addDistance(self):
        # dataframe中の 'wavelength', 'resolution_limit'を利用してカメラ長を計算する
        # 各数値は、self.df['wavelength'], self.df['resolution_limit']で取得できるが文字列の可能性があるので数値にしてから利用する
        self.df['wavelength'] = self.df['wavelength'].astype(float)
        self.df['resolution_limit'] = self.df['resolution_limit'].astype(float)
        self.df['dist_ds'] = self.df.apply(lambda x: self.calcDist(x['wavelength'], x['resolution_limit']), axis=1)
        # resolution limit は beamline.iniから読み込む
        # self.config : section=experiment, option=resol_raster
        roi_value = self.config.getint("experiment", "raster_roi", fallback=0)
        # roi flag
        # roi_value =1 -> roi_flag=True
        # roi_value =0 -> roi_flag=False
        resol_raster = self.config.getfloat("experiment", "resol_raster")
        if roi_value == 1:
            self.df['dist_raster'] = self.df.apply(lambda x: self.calcDist(x['wavelength'], resol_raster, True), axis=1)
        else:
            self.df['dist_raster'] = self.df.apply(lambda x: self.calcDist(x['wavelength'], resol_raster, False), axis=1)

    def makeCondList(self):
        # DataFrameとしてExcelファイルを読み込む　 →　self.df
        self.read_new()
        # self.dfの中にある 'puckid' の情報を展開する
        self.expandCompressedPinInfo()
        # 液体窒素ぶっ掛けの情報を管理してCSV用の情報へ変換
        self.checkLN2flag()
        # カメラのZoomに関する情報を管理してCSV用の情報へ変換
        self.checkZoomFlag()
        # カメラ長に関する情報をCSV用の情報へ変換
        self.addDistance()
        # ピンの温めに関する情報を管理してCSV用の情報へ変換
        self.checkPinFlag()
        self.splitBeamsizeInfo()
        # beam sizeから beamsize.config → Fluxを読み込んでDataFrameに入れる
        self.fillFlux()
        # default parameterをDataFrameに入れていく（beamline.iniからほとんど読み込んでいる）
        self.setDefaults()
        # raster scanの露光条件の決定
        self.defineScanCondition()
        # 露光条件について検討。transmission > 100% のときに露光時間とtransmissionを編集する
        self.modifyExposureConditions()
        # 結晶サイズについてのWarning（今は multi だけ)
        self.sizeWarning()

        # self.dfの内容をCSVファイルに書き出す
        # column の並び順は以下のように変更する
        # root_dir,p_index,mode,puckid,pinid,sample_name,wavelength,raster_vbeam,raster_hbeam,att_raster,
        # hebi_att,exp_raster,dist_raster,loopsize,score_min,score_max,maxhits,total_osc,osc_width,ds_vbeam,ds_hbeam,
        # exp_ds,dist_ds,dose_ds,offset_angle,reduced_fact,ntimes,meas_name,cry_min_size_um,cry_max_size_um,
        # hel_full_osc,hel_part_osc, raster_roi, ln2_flag, cover_scan_flag, zoomcap_flag, warm_time 
        # その他の値は self.df から読み込む
        self.columns = ['root_dir', 'p_index', 'mode', 'puckid', 'pinid', 'sample_name', 'wavelength', 'raster_vbeam', 'raster_hbeam', 'att_raster', \
                        'hebi_att', 'exp_raster', 'dist_raster', 'loopsize', 'score_min', 'score_max', 'maxhits', 'total_osc', 'osc_width', 'ds_vbeam', 'ds_hbeam', \
                        'exp_ds', 'dist_ds', 'dose_ds', 'offset_angle', 'reduced_fact', 'ntimes', 'meas_name', 'cry_min_size_um', 'cry_max_size_um', \
                        'hel_full_osc', 'hel_part_osc', 'raster_roi', 'ln2_flag', 'cover_scan_flag', 'zoomcap_flag', 'warm_time']       

        # ここで変数の型を明示的に指定する
        # float を想定しているもののみ
        # wavelength, resolution_limit, max_crystal_size, total_osc, osc_width
        # については float として読み込む
        set_types = {
            'wavelength': float,
            'raster_vbeam': float,
            'raster_hbeam': float,
            'att_raster': float,
            'hebi_att': float,
            'exp_raster': float,
            'dist_raster': float,
            'loopsize': float,
            'score_min': int,
            'score_max': int,
            'maxhits': int,
            'total_osc': float,
            'osc_width': float,
            'ds_vbeam': float,
            'ds_hbeam': float,
            'exp_ds': float,
            'dist_ds': float,
            'dose_ds': float,
            'offset_angle': float,
            'reduced_fact': float,
            'ntimes': int,
            'cry_min_size_um': float,
            'cry_max_size_um': float,
            'hel_full_osc': float,
            'hel_part_osc': float,
            'warm_time': float,
            'resolution_limit': float,
            'max_crystal_size': float,
            'total_osc': float,
            'osc_width': float,
        }
        # 型を指定する
        self.df = self.df.astype(set_types)

        # floatのフォーマットを指定
        float_format = '%.5f'
        # to_csv()メソッドでファイルに書き出す際にfloatのフォーマットを指定して書き出す
        zoo_csv_name = f"{self.csv_prefix}.csv"
        self.df.to_csv(zoo_csv_name, columns=self.columns, index=False, float_format=float_format)

        # 全パラメータの型を出力
        self.logger.info(f"Data types of all parameters in the DataFrame: {self.df.dtypes}")

if __name__ == "__main__":
    root_dir = os.getcwd()
    u2db = UserESA(sys.argv[1], root_dir, beamline="BL32XU")
    u2db.makeCondList()
    #u2db.read_new()
    #newdf = u2db.expandCompressedPinInfo()
    # CSV ファイルに書き出す
    #newdf.to_csv("check.csv", index=False)
    # u2db.df['ppf_raster']を 指数表記で出力
    #pd.options.display.float_format = '{:.2e}'.format
    #print(u2db.df['dist_raster'])
    #print(u2db.df['dist_ds'])
