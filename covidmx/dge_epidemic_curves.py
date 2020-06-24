from mapsmx import MapsMX
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import add
import os
import covidmx.utils_mult as utM
from numpy.core.arrayprint import _line_width
import covidmx.dge_plot as dgePlt


class DGEEpidemicCurves(dgePlt.DGEBase):
    """
    Class to plot weekly multipliers dge information
    """

    def plotStringencyDates(self, weekly=True):
        susy=self.stringency_dates["susana"]-self.dias_confirmados[self.discardFirstDays]
        nueva_norm=self.stringency_dates["nueva_norm"]-self.dias_confirmados[self.discardFirstDays]        
        
        if weekly:
            plt.axvline(x=susy.days/7,  linewidth=4, c=(1.0, 0, 0,0.5))
            plt.axvline(x=nueva_norm.days/7,  linewidth=4, c=(1, 0.65, 0,0.5))
        else:
            plt.axvline(x=susy.days,  linewidth=4, c=(1.0, 0, 0,0.5))
            plt.axvline(x=nueva_norm.days,  linewidth=4, c=(1, 0.65, 0,0.5))    
        return susy, nueva_norm

    def plotHistoricDaily(self, fig, per_day, metric, historic_days, plot_position, line_color='bo-', line_width=0.4, log_scale=True ):        
        xAxis=[x for x in range(len(historic_days))]
        ax = fig.add_subplot(plot_position);
        _=plt.xticks(xAxis, historic_days , rotation='vertical');ax.grid('on')#plt.grid();#plt.show()#[x[0][8:13] for x in historic_days]
        if log_scale: ax.set_yscale('log'); 
        ax.yaxis.grid(b=True, which='both', linestyle='--')
        if metric == '':lMetric=""
        else: lMetric="{} por dìa".format(metric)
        ax.plot( xAxis,per_day, line_color, label=lMetric, linewidth=line_width, markersize=line_width*4.5 );
        ax.set_ylabel(r"nùmero {}".format(metric)); ax.legend(fontsize="small") #, loc=6
        
        susy, _= self.plotStringencyDates( weekly=False)
        for ith_monday in range((susy.days)%7,len(xAxis),7):
            plt.axvline(x=ith_monday, linewidth=1,linestyle='dashed')
        uncertain_days=self.delayReportDays-self.discardLatestDays
        if uncertain_days>0:
            plt.axvline(x=len(historic_days)-uncertain_days,  linewidth=4, c=(1.0, 0, 0,0.5), linestyle='dashed')
            
    def plotHistoricWeekly(self, fig, muertos_week, multipliers, muertos_week_suspect, muertos_mult_suspect, metric, historic_week, plot_position ):
        ax = fig.add_subplot(220+plot_position); xAxis=[x for x in range(len(historic_week))]
        _=plt.xticks(xAxis, historic_week , rotation='vertical');ax.grid('on');#plt.show()#[x[0][8:13] for x in historic_week]
        ax.yaxis.grid(b=True, which='minor', linestyle='--'); #ax.set_yscale('log'); 
        ax.plot(xAxis,muertos_week,'o-', label="{} per week".format(metric) );
        
        ax.plot(xAxis,muertos_week_suspect,'mo-', label="{} + sospechosos per week".format(metric) );
        
        ax.set_ylabel(r"{} number".format(metric)); ax.legend(fontsize="small")#, loc=6
        self.plotStringencyDates()
        
        ax2 = fig.add_subplot(222+plot_position);#ax2 = ax.twinx(); 
        _=plt.xticks(xAxis, historic_week , rotation='vertical');ax2.grid('on');
        ax2.plot(xAxis,multipliers,'^-', label="Multipliers")
        ax2.plot(xAxis,muertos_mult_suspect,'m^-', label="Multipliers + sospechosos")
        ax2.yaxis.grid(b=True, which='minor', linestyle='--'); ax2.set_ylim((0,2.0))# plt.ylim(4.0)
        ax2.set_ylabel(r"{} multiplier".format(metric)); ax2.legend(fontsize="small");ax.grid('on') #, loc=5
        
#         susy=self.stringency_dates["susana"]-self.dias_confirmados[self.discardFirstDays]
#         nueva_norm=self.stringency_dates["nueva_norm"]-self.dias_confirmados[self.discardFirstDays]    
        self.plotStringencyDates()
        uncertain_days=self.delayReportDays-self.discardLatestDays
        if uncertain_days>0:
            ax.axvline(x=len(historic_week)-uncertain_days/7.0,  linewidth=4, c=(1.0, 0, 0,0.5), linestyle='dashed')
            ax2.axvline(x=len(historic_week)-uncertain_days/7.0,  linewidth=4, c=(1.0, 0, 0,0.5), linestyle='dashed')

    def casesPerDayInRange(self, muertos_df, casos_df, muertos_sospechosos_df, sospechosos_df):
        casos_dia=[]; muertos_dia=[]; casos_dia_ingreso=[]; historic_days=[];historic_days_full=[]#  #muertos_por_dia=[];  
        casos_d_s=[]; muertos_d_s=[]; casos_d_s_ingreso=[]
        upper_bound= len(self.dias_confirmados)-1-self.discardLatestDays#TODO: compute 7 days avg; get data/numbers from -/+3 days # date_format='%m-%d'; 
        for day in pd.date_range(self.dias_confirmados[self.discardFirstDays],self.dias_confirmados[upper_bound]): #dias_muertes  
            muertos_dia, casos_dia, casos_dia_ingreso=utM.casosPorDia(muertos_df, casos_df, day, muertos_dia, casos_dia, casos_dia_ingreso)
            muertos_d_s, casos_d_s, casos_d_s_ingreso=utM.casosPorDia(muertos_sospechosos_df, sospechosos_df, day, muertos_d_s, casos_d_s, casos_d_s_ingreso)
            if len(historic_days)==0 or day.date().day==1:
                historic_days.append( '{}/{}'.format(day.date().day,day.date().month) )
            else:
                historic_days.append( '{}'.format(day.date().day) )
            historic_days_full.append( '{}/{}'.format(day.date().day,day.date().month) )
        
        return  casos_dia, muertos_dia, casos_dia_ingreso, [historic_days,historic_days_full], casos_d_s, muertos_d_s, casos_d_s_ingreso
    
    

    def plotHistoric(self, state=None, municipality=None, ploat_all=False):
        if state is not None: 
            plot_data = self.dge_data[ self.dge_data['entidad_res'].str.lower() == state.lower() ] #plot_data['municipio_res'].unique()
            
            population=self.muniDF[ self.muniDF["NOM_ENT"].str.lower()==state.lower() ]["POB_TOT"].sum()
            if municipality is not None:
#                 plot_data['municipio_res'].unique() #TODO: by pass bad namening find more similar municipality.lower(), verify length ofunique municipalities
                plot_data = plot_data[ plot_data['municipio_res'].str.lower() == municipality.lower() ] 
                state=municipality+", "+state #plot_data.keys()# 
                
                cve_mun=plot_data["cve_mun"].unique()[0].split('_')
                if len(cve_mun[1])<3: cve_mun[1]='0'+cve_mun[1]
                if len(cve_mun[1])<3: cve_mun[1]='0'+cve_mun[1]
                cve_mun=int( cve_mun[0]+cve_mun[1] )
                population=self.muniDF[ self.muniDF["CVE_MUN"]==cve_mun ]["POB_TOT"].values[0]
            
        else:
            plot_data = self.dge_data
            state="Nacional"
            population=self.muniDF["POB_TOT"].sum()
            
        metric="muertos"; metricCasos="casos"
        casos_df=plot_data[ plot_data['resultado']=='confirmados' ];         
        muertos_df = casos_df[casos_df['muertos']==1]#casos_df['muertos'].sum()        
        self.dias_confirmados=np.sort( casos_df['fecha_sintomas'].unique() );#dias_muertes=np.sort( muertos_df['fecha_def'].unique() )  
        # np.sort( casos_df['fecha_ingreso'].unique() )                  
        
        sospechosos_df=plot_data[ plot_data['resultado']=='sospechosos' ]
        muertos_sospechosos_df = sospechosos_df[sospechosos_df['muertos']==1]#casos_df['muertos'].sum()        
        self.dias_sospechosos=np.sort( sospechosos_df['fecha_sintomas'].unique() );
        
        self.discardLatestDays=0;self.discardFirstDays=0; self.delayReportDays=14
        casos_dia, muertos_dia, casos_dia_ingreso, historic_d, casos_d_s, muertos_d_s, casos_d_s_ingreso = self.casesPerDayInRange( muertos_df, casos_df, muertos_sospechosos_df, sospechosos_df)
        if ploat_all:
#             title="Muertes y casos (fecha de sìntomas/ingreso) daily historic" ;fig = plt.figure(figsize=(20.0, 15.0));            
#             fig.suptitle("{}\n {}".format(title, state) )
    #         fig.savefig(os.path.join(mobiVisuRes,title.replace(' ', '_')+".png"), bbox_inches='tight')            
            historic_days=historic_d[0]
            title="Muertes y casos (confirmados, sospechosos y con+sos fecha de sìntomas/ingreso) daily historic" ; fig = plt.figure(figsize=(20.0, 15.0));
            self.plotHistoricDaily(fig, muertos_dia, metric, historic_days, 311 ) #TODO: Factorize; verify symptoms start
            self.plotHistoricDaily(fig, list( map(add, muertos_dia, muertos_d_s) ), metric+" con+sos", historic_days, 311, line_color='mo-' )
            self.plotHistoricDaily(fig, muertos_d_s, metric, historic_days, 311, line_color='ro-', line_width=0.6)            
            
            self.plotHistoricDaily(fig, casos_dia, metricCasos, historic_days, 312 )
            self.plotHistoricDaily(fig, list( map(add, casos_dia, casos_d_s) ), metricCasos+" con+sos", historic_days, 312, line_color='mo-'  )   
            self.plotHistoricDaily(fig, casos_d_s, metricCasos+" sos", historic_days, 312, line_color='ro-', line_width=0.6 )
            
            self.plotHistoricDaily(fig, casos_dia_ingreso, metricCasos+" ingreso sos", historic_days, 313 )       
            self.plotHistoricDaily(fig, list( map(add, casos_dia_ingreso, casos_d_s_ingreso) ), metricCasos+" ingreso con+sos", historic_days, 313, line_color='mo-' ) #TODO: verify symptoms start
            self.plotHistoricDaily(fig, casos_d_s_ingreso, metricCasos+" ingreso", historic_days, 313, line_color='ro-', line_width=0.6 ) 
            fig.suptitle("{}\n {}".format(title, state) )
    #         fig.savefig(os.path.join(mobiVisuRes,title.replace(' ', '_')+".png"), bbox_inches='tight')
        
        self.discardFirstDays=17
        casos_dia, muertos_dia, casos_dia_ingreso, historic_d, casos_d_s, muertos_d_s, casos_d_s_ingreso = self.casesPerDayInRange( muertos_df, casos_df, muertos_sospechosos_df, sospechosos_df)
#         casos_dia=[]; muertos_dia=[]; casos_dia_ingreso=[]; historic_days=[]; casos_d_s=[]; muertos_d_s=[]; casos_d_s_ingreso=[]
#         upper_bound= len(self.dias_confirmados)-1-self.discardLatestDays#TODO: compute 7 days avg; get data/numbers from -/+3 days # date_format='%m-%d'; 
#         for day in pd.date_range(self.dias_confirmados[self.discardFirstDays],self.dias_confirmados[upper_bound]): #dias_muertes  
#             muertos_dia, casos_dia, casos_dia_ingreso=utM.casosPorDia(muertos_df, casos_df, day, muertos_dia, casos_dia, casos_dia_ingreso)
#             muertos_d_s, casos_d_s, casos_d_s_ingreso=utM.casosPorDia(muertos_sospechosos_df, sospechosos_df, day, muertos_d_s, casos_d_s, casos_d_s_ingreso)
#             historic_days.append( '{}/{}'.format(day.date().day,day.date().month) )
        
        muertos_plus_suspect = list( map(add, muertos_dia, muertos_d_s) )
        casos_plus_suspect=list( map(add, casos_dia, casos_d_s) )
        if ploat_all:
            historic_days=historic_d[0]
            title="Muertes (confirmados+sospechosos) daily historic"; fig = plt.figure(figsize=(20.0, 15.0));  
            muertos_dia_ma = utM.smooth(np.array(muertos_dia), window_len=7).tolist()[3:-3]; metricSmooth="" # 
#             muertos_dia_ma2 = utM.smooth(np.array(muertos_dia), window_len=7).tolist()
#             muertos_dia_ma = muertos_dia[:3]+np.convolve(muertos_dia, np.ones((7,))/7, mode='valid').tolist()+muertos_dia[-3
            self.plotHistoricDaily(fig, muertos_dia, metric, historic_days, 111, log_scale=False )            
            self.plotHistoricDaily(fig, muertos_dia_ma, metricSmooth, historic_days, 111, line_color='b-',log_scale=False , line_width=1.2)

            muertos_plus_suspect_ma = utM.smooth(np.array(muertos_plus_suspect), window_len=7).tolist()[3:-3] #            
            self.plotHistoricDaily(fig, muertos_plus_suspect, metric+" (confirmados+sospechosos)", historic_days, 111, line_color='mo-', log_scale=False)
            self.plotHistoricDaily(fig, muertos_plus_suspect_ma, metricSmooth, historic_days, 111, line_color='m-',log_scale=False, line_width=1.2 )
            fig.suptitle("{}\n {}".format(title, state) )            
            
            title="Casos (confirmados+sospechosos fecha de sìntomas/ingreso) daily historic"; fig = plt.figure(figsize=(20.0, 15.0));  
            
            casos_dia_ma = utM.smooth(np.array(casos_dia), window_len=7, window='flat').tolist()[3:-3] #
            self.plotHistoricDaily(fig, casos_dia, metricCasos, historic_days, 111, log_scale=False )
            self.plotHistoricDaily(fig, casos_dia_ma, metricSmooth, historic_days, 111, line_color='b-',log_scale=False, line_width=1.2 )
            
            casos_plus_suspect_ma = utM.smooth(np.array(casos_plus_suspect), window_len=7, window='flat').tolist()[3:-3] #            
            self.plotHistoricDaily(fig, casos_plus_suspect, metricCasos+" (confirmados+sospechosos)", historic_days, 111, line_color='mo-', log_scale=False)        
            self.plotHistoricDaily(fig, casos_plus_suspect_ma, metricSmooth, historic_days, 111, line_color='m-',log_scale=False, line_width=1.2 )
            
    #         self.plotHistoricDaily(fig, list( map(add, casos_dia_ingreso, casos_d_s_ingreso) ), metricCasos+" ingreso", historic_days, 313 ) #TODO: verify symptoms start
    #         self.plotHistoricDaily(fig, casos_d_s_ingreso, metricCasos+" ingreso", historic_days, 313, line_color='ro-' )
            fig.suptitle("{}\n {}".format(title, state) )
        
        historic_week=[]; muertos_week=[];  muertos_week_mult=[0]; casos_week=[];  casos_week_mult=[0]; #TODO: weekly from latest day 
        historic_days=historic_d[1]; muertos_week_suspect=[];  muertos_mult_suspect=[0]; casos_week_suspect=[];  casos_mult_suspect=[0];
        for wk in range( int((len(muertos_dia) )/7) ):           #TODO: verify metric in FT.com
            muertos_week.append( np.sum( muertos_dia[wk*7:wk*7+7] ) ); #muertos_dia[wk*7:wk*7+7];wk*7;wk*7+7
            casos_week.append( np.sum( casos_dia[wk*7:wk*7+7] ) ); 
            
            muertos_week_suspect.append( np.sum( muertos_plus_suspect[wk*7:wk*7+7] ) ); #muertos_dia[wk*7:wk*7+7];wk*7;wk*7+7
            casos_week_suspect.append( np.sum( casos_plus_suspect[wk*7:wk*7+7] ) ); 
            
            historic_week.append( historic_days[wk*7] )
            if wk>0:                
                muertos_week_mult.append(muertos_week[-1]/muertos_week[-2])
                casos_week_mult.append(casos_week[-1]/casos_week[-2])  
                
                muertos_mult_suspect.append(muertos_week_suspect[-1]/muertos_week_suspect[-2])
                casos_mult_suspect.append(casos_week_suspect[-1]/casos_week_suspect[-2])

        muertos_week_mult=utM.cleanMult(muertos_week_mult); casos_week_mult=utM.cleanMult(casos_week_mult)    
        
        casosT=np.sum(casos_week_suspect); muertosT=np.sum(muertos_week_suspect);casosAct=np.sum(casos_week[-3:]); casosActSusp=np.sum(casos_week_suspect[-3:])
        casosActInci=casosAct*100000.0/population; casosActSuspInci=casosActSusp*100000.0/population
        title="Casos {} ({}) y muertos {} ({}) (confirmados+sospechosos)".format(np.sum(casos_week), casosT, np.sum(muertos_week), muertosT ); 
        title2="\nCasos activos last 3weeks {} ({}) Incidencia {:.2f} ({:.2f}) casos/100k ".format(casosAct, casosActSusp, casosActInci, casosActSuspInci); 
        fig = plt.figure(figsize=(10.0, 10.0));                
        self.plotHistoricWeekly( fig, casos_week, casos_week_mult, casos_week_suspect, casos_mult_suspect, metricCasos, historic_week, 1 )        
        self.plotHistoricWeekly( fig, muertos_week, muertos_week_mult, muertos_week_suspect, muertos_mult_suspect, metric, historic_week, 2 )
        casosInci=casosT*100000.0/population; muertosInci=muertosT*100000.0/population
        fig.suptitle("{}\nIncidencia {:.2f} casos/100k  y {:.2f} muertos/100k  Pob{} {}\n{} weekly and multipliers".format(title, casosInci, muertosInci, population, title2, state) )

#         fig.suptitle("{}\n {}".format(title, metricRangeStr) )
        fig.savefig(os.path.join("/data/covid/fer",state.replace(", ","_")+".png"), bbox_inches='tight')

