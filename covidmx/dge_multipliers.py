from mapsmx import MapsMX
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import add

import covidmx.utils_mult as utM

class DGEMultipliers:
    """
    Class to plot weekly multipliers dge information
    """

    def __init__(self, dge_data, catalogue, description):

        self.dge_data = self.prepare_data(dge_data)
        self.dge_data['cve_ent'] = self.dge_data['cve_ent'].astype(int).astype(str)
        self.catalogue = catalogue
        self.description = description

        #Downloading geo information
        state_geo = MapsMX().get_geo('state')
        state_geo['cve_ent'] = state_geo['cve_ent'].astype(int).astype(str)

        mun_geo = MapsMX().get_geo('municipality')
        mun_geo[['cve_ent', 'cve_mun']] = mun_geo[['cve_ent', 'cve_mun']].astype(int).astype(str)
        mun_geo['cve_mun'] = mun_geo['cve_ent'] + '_' + mun_geo['cve_mun']


        self.state_geo = state_geo
        self.mun_geo = mun_geo
        self.available_states = self.dge_data['entidad_res'].unique()
        self.available_status = ['confirmados', 'negativos', 'sospechosos', 'muertos']



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

    def plotHistoricDaily(self, fig, per_day, metric, historic_days, plot_position, line_color='bo-', log_scale=True ):        
        xAxis=[x for x in range(len(historic_days))]
        ax = fig.add_subplot(plot_position);
        _=plt.xticks(xAxis, historic_days , rotation='vertical');ax.grid('on')#plt.grid();#plt.show()#[x[0][8:13] for x in historic_days]
        if log_scale: ax.set_yscale('log'); 
        ax.yaxis.grid(b=True, which='both', linestyle='--')
        ax.plot(xAxis,per_day, line_color, label="{} por dìa".format(metric) );
        ax.set_ylabel(r"{} nùmero".format(metric)); ax.legend(fontsize="small", loc=6)
        
        susy, _= self.plotStringencyDates( weekly=False)
        for ith_monday in range((susy.days)%7,len(xAxis),7):
            plt.axvline(x=ith_monday, linewidth=1,linestyle='dashed')
        uncertain_days=self.delayReportDays-self.discardLatestDays
        if uncertain_days>0:
            plt.axvline(x=len(historic_days)-uncertain_days,  linewidth=4, c=(1.0, 0, 0,0.5), linestyle='dashed')
            
    def plotHistoricWeekly(self, fig, muertos_week, multipliers, metric, historic_week, plot_position ):
        ax = fig.add_subplot(220+plot_position); xAxis=[x for x in range(len(historic_week))]
        _=plt.xticks(xAxis, historic_week , rotation='vertical');ax.grid('on');#plt.show()#[x[0][8:13] for x in historic_week]
        ax.yaxis.grid(b=True, which='minor', linestyle='--'); #ax.set_yscale('log'); 
        ax.plot(xAxis,muertos_week,'o-', label="{} per week".format(metric) );
        ax.set_ylabel(r"{} number".format(metric)); ax.legend(fontsize="small", loc=6)
        self.plotStringencyDates()
        
        ax2 = fig.add_subplot(222+plot_position);#ax2 = ax.twinx(); 
        _=plt.xticks(xAxis, historic_week , rotation='vertical');ax2.grid('on');
        ax2.plot(xAxis,multipliers,'o^-', label="Multipliers")
        ax2.yaxis.grid(b=True, which='minor', linestyle='--'); ax2.set_ylim((0,2.0))# plt.ylim(4.0)
        ax2.set_ylabel(r"{} multiplier".format(metric)); ax2.legend(fontsize="small", loc=5);ax.grid('on')
        
#         susy=self.stringency_dates["susana"]-self.dias_confirmados[self.discardFirstDays]
#         nueva_norm=self.stringency_dates["nueva_norm"]-self.dias_confirmados[self.discardFirstDays]    
        self.plotStringencyDates()
        uncertain_days=self.delayReportDays-self.discardLatestDays
        if uncertain_days>0:
            ax.axvline(x=len(historic_week)-uncertain_days/7.0,  linewidth=4, c=(1.0, 0, 0,0.5), linestyle='dashed')
            ax2.axvline(x=len(historic_week)-uncertain_days/7.0,  linewidth=4, c=(1.0, 0, 0,0.5), linestyle='dashed')

    def casesPerDayInRange(self, muertos_df, casos_df, muertos_sospechosos_df, sospechosos_df):
        casos_dia=[]; muertos_dia=[]; casos_dia_ingreso=[]; historic_days=[];#  #muertos_por_dia=[];  
        casos_d_s=[]; muertos_d_s=[]; casos_d_s_ingreso=[]
        upper_bound= len(self.dias_confirmados)-1-self.discardLatestDays#TODO: compute 7 days avg; get data/numbers from -/+3 days # date_format='%m-%d'; 
        for day in pd.date_range(self.dias_confirmados[self.discardFirstDays],self.dias_confirmados[upper_bound]): #dias_muertes  
            muertos_dia, casos_dia, casos_dia_ingreso=utM.casosPorDia(muertos_df, casos_df, day, muertos_dia, casos_dia, casos_dia_ingreso)
            muertos_d_s, casos_d_s, casos_d_s_ingreso=utM.casosPorDia(muertos_sospechosos_df, sospechosos_df, day, muertos_d_s, casos_d_s, casos_d_s_ingreso)
            historic_days.append( '{}/{}'.format(day.date().day,day.date().month) )
        
        return  casos_dia, muertos_dia, casos_dia_ingreso, historic_days, casos_d_s, muertos_d_s, casos_d_s_ingreso
    
    

    def plotHistoric(self, state=None):
        if state is not None: 
            plot_data = self.dge_data[self.dge_data['entidad_res'].str.lower() == state.lower()]
        else:
            plot_data = self.dge_data
            
        metric="muertos"; metricCasos="casos"
        casos_df=plot_data[ plot_data['resultado']=='confirmados' ];         
        muertos_df = casos_df[casos_df['muertos']==1]#casos_df['muertos'].sum()        
        self.dias_confirmados=np.sort( casos_df['fecha_sintomas'].unique() );#dias_muertes=np.sort( muertos_df['fecha_def'].unique() )  
        # np.sort( casos_df['fecha_ingreso'].unique() )                  
        
        sospechosos_df=plot_data[ plot_data['resultado']=='sospechosos' ]
        muertos_sospechosos_df = sospechosos_df[sospechosos_df['muertos']==1]#casos_df['muertos'].sum()        
        self.dias_sospechosos=np.sort( sospechosos_df['fecha_sintomas'].unique() );
        
        self.discardLatestDays=0;self.discardFirstDays=0; self.delayReportDays=14
        casos_dia, muertos_dia, casos_dia_ingreso, historic_days, casos_d_s, muertos_d_s, casos_d_s_ingreso = self.casesPerDayInRange( muertos_df, casos_df, muertos_sospechosos_df, sospechosos_df)

        title="Muertes y casos (fecha de sìntomas/ingreso) daily historic" ;fig = plt.figure(figsize=(20.0, 15.0));
        self.plotHistoricDaily(fig, muertos_dia, metric, historic_days, 311 )
        self.plotHistoricDaily(fig, casos_dia, metricCasos, historic_days, 312 )
        self.plotHistoricDaily(fig, casos_dia_ingreso, metricCasos+" ingreso", historic_days, 313 ) #TODO: verify symptoms start
        fig.suptitle("{}\n {}".format(title, state) )
#         fig.savefig(os.path.join(mobiVisuRes,title.replace(' ', '_')+".png"), bbox_inches='tight')            
        
        title="Muertes y casos (confirmados+sospechosos fecha de sìntomas/ingreso) daily historic" ; fig = plt.figure(figsize=(20.0, 15.0));
        self.plotHistoricDaily(fig, list( map(add, muertos_dia, muertos_d_s) ), metric, historic_days, 311 )
        self.plotHistoricDaily(fig, muertos_d_s, metric, historic_days, 311, line_color='ro-')        
        self.plotHistoricDaily(fig, list( map(add, casos_dia, casos_d_s) ), metricCasos, historic_days, 312 )
        self.plotHistoricDaily(fig, casos_d_s, metricCasos, historic_days, 312, line_color='ro-' )        
        self.plotHistoricDaily(fig, list( map(add, casos_dia_ingreso, casos_d_s_ingreso) ), metricCasos+" ingreso", historic_days, 313 ) #TODO: verify symptoms start
        self.plotHistoricDaily(fig, casos_d_s_ingreso, metricCasos+" ingreso", historic_days, 313, line_color='ro-' )        
        fig.suptitle("{}\n {}".format(title, state) )
#         fig.savefig(os.path.join(mobiVisuRes,title.replace(' ', '_')+".png"), bbox_inches='tight')
        
        self.discardFirstDays=17
        casos_dia, muertos_dia, casos_dia_ingreso, historic_days, casos_d_s, muertos_d_s, casos_d_s_ingreso = self.casesPerDayInRange( muertos_df, casos_df, muertos_sospechosos_df, sospechosos_df)
#         casos_dia=[]; muertos_dia=[]; casos_dia_ingreso=[]; historic_days=[]; casos_d_s=[]; muertos_d_s=[]; casos_d_s_ingreso=[]
#         upper_bound= len(self.dias_confirmados)-1-self.discardLatestDays#TODO: compute 7 days avg; get data/numbers from -/+3 days # date_format='%m-%d'; 
#         for day in pd.date_range(self.dias_confirmados[self.discardFirstDays],self.dias_confirmados[upper_bound]): #dias_muertes  
#             muertos_dia, casos_dia, casos_dia_ingreso=utM.casosPorDia(muertos_df, casos_df, day, muertos_dia, casos_dia, casos_dia_ingreso)
#             muertos_d_s, casos_d_s, casos_d_s_ingreso=utM.casosPorDia(muertos_sospechosos_df, sospechosos_df, day, muertos_d_s, casos_d_s, casos_d_s_ingreso)
#             historic_days.append( '{}/{}'.format(day.date().day,day.date().month) )

        title="Muertes (confirmados+sospechosos fecha de sìntomas/ingreso) daily historic"; fig = plt.figure(figsize=(20.0, 15.0));  
        self.plotHistoricDaily(fig, muertos_dia, metric, historic_days, 111, log_scale=False )
        self.plotHistoricDaily(fig, list( map(add, muertos_dia, muertos_d_s) ), metric, historic_days, 111, line_color='ro-', log_scale=False)
        fig.suptitle("{}\n {}".format(title, state) )
        
        title="Casos (confirmados+sospechosos fecha de sìntomas/ingreso) daily historic"; fig = plt.figure(figsize=(20.0, 15.0));  
        self.plotHistoricDaily(fig, casos_dia, metricCasos, historic_days, 111, log_scale=False )
        self.plotHistoricDaily(fig, list( map(add, casos_dia, casos_d_s) ), metricCasos, historic_days, 111, line_color='ro-', log_scale=False)        
#         self.plotHistoricDaily(fig, list( map(add, casos_dia_ingreso, casos_d_s_ingreso) ), metricCasos+" ingreso", historic_days, 313 ) #TODO: verify symptoms start
#         self.plotHistoricDaily(fig, casos_d_s_ingreso, metricCasos+" ingreso", historic_days, 313, line_color='ro-' )
        fig.suptitle("{}\n {}".format(title, state) )
        
        
        historic_week=[]; muertos_week=[];  muertos_week_mult=[0]; casos_week=[];  casos_week_mult=[0];
        muertos_week_suspect=[];  muertos_mult_suspect=[0]; casos_week_suspect=[];  casos_mult_suspect=[0];
        for wk in range( int((len(muertos_dia) )/7) ):           #TODO: verify metric in FT.com
            muertos_week.append( np.sum( muertos_dia[wk*7:wk*7+7] ) ); #muertos_dia[wk*7:wk*7+7];wk*7;wk*7+7
            casos_week.append( np.sum( casos_dia[wk*7:wk*7+7] ) ); 
            historic_week.append( historic_days[wk*7] )
            if wk>0:                
                muertos_week_mult.append(muertos_week[-1]/muertos_week[-2])
                casos_week_mult.append(casos_week[-1]/casos_week[-2])
        
        muertos_week_mult=utM.cleanMult(muertos_week_mult); casos_week_mult=utM.cleanMult(casos_week_mult)    

        title="Casos (confirmados+sospechosos fecha de sìntomas/ingreso) weekly and multipliers historic"; fig = plt.figure(figsize=(10.0, 10.0));                
        self.plotHistoricWeekly( fig, muertos_week, muertos_week_mult.tolist(), metric, historic_week, 2 )
        self.plotHistoricWeekly( fig, casos_week, casos_week_mult.tolist(), metricCasos, historic_week, 1 )
        
#         fig.suptitle("{}\n {}".format(title, metricRangeStr) )
#         fig.savefig(os.path.join(mobiVisuRes,title.replace(' ', '_')+".png"), bbox_inches='tight')
        
        

    def prepare_data(self, df):

        df = df.rename(columns={
                       'entidad_res_original': 'cve_ent',
                       'municipio_res_original': 'cve_mun'
                       })
        
        df['muertos'] = df['fecha_def'].notna().astype(int)

        replace_resultado = {'Positivo SARS-CoV-2': 'confirmados',
                             'No positivo SARS-CoV-2': 'negativos',
                             'Resultado pendiente':'sospechosos'}

        df['resultado'] = df['resultado'].replace(replace_resultado)
        df = pd.concat([df, pd.get_dummies(df['resultado'])], axis=1)

        int_vars = list(replace_resultado.values()) + ['muertos']
        df[int_vars] = df[int_vars].astype(int) #19764 09-07 
        
#         self.plotHistoric(df)

        return df

    def plot_map(self, status='confirmados', state=None,
                 add_municipalities=False, save_file_name = None,
                 cmap='Reds',
                 scheme='quantiles', k=4, legend=True, zorder=1,
                 missing_kwds={'color': 'lightgray', 'label': 'Sin info'}, **kwargs):
        """
        Plot geography information

        Parameters
        ----------
        status: str
            One of confirmados, sospechosos, negativos, muertos
        state: str
            Plot particular state.
        add_municipalities: bool
            Wheter add municipalities to plot
        """

        assert status in self.available_status, 'Please provide some of the following status: {}'.format(', '.join(self.available_status))
        if state is not None:
            assert state in self.available_states, 'Please provide some of the following states: {}'.format(', '.join(self.available_states))

        # if last_date_to_consider is not None:
        #     last_date = pd.to_datetime(last_date_to_consider, format=format_date)
        #     plot_data = self.dge_data[self.dge_data['fecha_sintomas']<=last_date]
        # else:
        #     plot_data = self.dge_data

        group_cols = ['entidad_res', 'cve_ent']

        if add_municipalities:
            group_cols += ['municipio_res', 'cve_mun']

        needed_cols = [status] + group_cols

        plot_data = self.dge_data[needed_cols]
        state_geo_plot = self.state_geo
        mun_geo_plot = self.mun_geo

        if state is not None:
            plot_data = plot_data[plot_data['entidad_res'].str.lower() == state.lower()]
            cve_ent = str(plot_data['cve_ent'].unique()[0])
            state_geo_plot = self.state_geo[self.state_geo['cve_ent']==cve_ent]
            mun_geo_plot = self.mun_geo[self.mun_geo['cve_ent']==cve_ent]

        plot_data = plot_data.groupby(group_cols).agg(sum).reset_index()

        if add_municipalities:
            plot_data = plot_data.drop(columns='cve_ent')
            plot_data = mun_geo_plot.merge(plot_data, how='left', on='cve_mun')

            geometry = 'geometry_mun'
        else:
            plot_data = state_geo_plot.merge(plot_data, how='left', on='cve_ent')
            geometry = 'geometry_ent'

        base = state_geo_plot.boundary.plot(color=None,
                                            edgecolor='black',
                                            linewidth=0.6,
                                            figsize=(10,9))

        if add_municipalities and state is not None:
            mun_geo_plot.boundary.plot(ax=base, color=None,
                                       edgecolor='black',
                                       linewidth=0.2)

        plot_obj = plot_data.set_geometry(geometry).plot(ax=base,
                                                         column=status,
                                                         cmap=cmap,
                                                         scheme=scheme,
                                                         k=k,
                                                         legend=legend,
                                                         zorder=zorder,
                                                         missing_kwds=missing_kwds,
                                                         **kwargs)
        base.set_axis_off()
        plt.axis('equal')

        title = 'Casos ' + status + ' por COVID-19'

        if state is not None:
            title += '\n'
            title += state.title()

        act_date = self.date
        if self.date is None:
            act_date = self.dge_data['fecha_actualizacion'][0]
        act_date = act_date.date()
        act_date = str(act_date)

        title += '\n'
        title += 'Fecha de actualizacion de los datos: {}'.format(act_date)


        plt.title(title, fontsize=20)


        if save_file_name is not None:
            plt.savefig(save_file_name, bbox_inches='tight', pad_inches=0)
            plt.close()
        else:
            plt.show()





        return plot_obj
