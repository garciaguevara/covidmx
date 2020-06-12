from mapsmx import MapsMX
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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





    def plotHistoricDaily(self, fig, num_per_day, metric, historic_days, plot_position ):        
        xAxis=[x for x in range(len(historic_days))]
        ax = fig.add_subplot(plot_position);
        _=plt.xticks(xAxis, historic_days , rotation='vertical');plt.grid();#plt.show()#[x[0][8:13] for x in historic_days]
        ax.set_yscale('log'); ax.yaxis.grid(b=True, which='minor', linestyle='--')
        ax.plot(xAxis,num_per_day,'o-', label="{} por dìa".format(metric) );
        ax.set_ylabel(r"{} nùmero".format(metric)); ax.legend(fontsize="small", loc=6)
        
        
    def plotHistoricWeekly(self, fig, muertos_week, multipliers, metric, historic_week, plot_position ):
        ax = fig.add_subplot(220+plot_position); xAxis=[x for x in range(len(historic_week))]
        _=plt.xticks(xAxis, historic_week , rotation='vertical');plt.grid();#plt.show()#[x[0][8:13] for x in historic_week]
        ax.yaxis.grid(b=True, which='minor', linestyle='--'); #ax.set_yscale('log'); 
        ax.plot(xAxis,muertos_week,'o-', label="{} per week".format(metric) );
        ax.set_ylabel(r"{} number".format(metric)); ax.legend(fontsize="small", loc=6)
        
        ax2 = fig.add_subplot(222+plot_position);#ax2 = ax.twinx(); 
        _=plt.xticks(xAxis, historic_week , rotation='vertical');
        ax2.plot(xAxis,multipliers,'r^-', label="Multipliers")
        ax2.yaxis.grid(b=True, which='minor', linestyle='--'); ax2.set_ylim((0,4.0))# plt.ylim(4.0)
        ax2.set_ylabel(r"{} multiplier".format(metric)); ax2.legend(fontsize="small", loc=5);plt.grid();plt.show()


    def prepare_data(self, df):

        df = df.rename(columns={
                       'entidad_res_original': 'cve_ent',
                       'municipio_res_original': 'cve_mun'
                       })
        metric="muertos"; metricConfirmados="confirmados"
        df['muertos'] = df['fecha_def'].notna().astype(int)

        replace_resultado = {'Positivo SARS-CoV-2': 'confirmados',
                             'No positivo SARS-CoV-2': 'negativos',
                             'Resultado pendiente':'sospechosos'}

        df['resultado'] = df['resultado'].replace(replace_resultado)
        df = pd.concat([df, pd.get_dummies(df['resultado'])], axis=1)

        int_vars = list(replace_resultado.values()) + ['muertos']
        df[int_vars] = df[int_vars].astype(int) #19764 09-07 
        confirmados_df=df[ df['resultado']=='confirmados' ];         
        muertos_df = confirmados_df[confirmados_df['muertos']==1]#confirmados_df['muertos'].sum()        
        dias_confirmados=np.sort( confirmados_df['fecha_sintomas'].unique() );#dias_muertes=np.sort( muertos_df['fecha_def'].unique() )  
        # np.sort( confirmados_df['fecha_ingreso'].unique() )                  
        
        num_confirmados_dia=[]; num_muertos_dia=[]; historic_days=[]; discardLatestDays=9; num_confirmados_por_dia_ingreso=[] #muertos_por_dia=[];  
        date_format='%m-%d';  #TODO: compute 7 days avg; get data/numbers from -/+3 days
        for day in pd.date_range(dias_confirmados[0],dias_confirmados[-discardLatestDays]): #dias_muertes  
            muertos_por_dia = muertos_df[ muertos_df['fecha_def']==day.date().__str__() ] #.append(
            num_muertos_dia.append( len(muertos_por_dia) ) #[-1]            
            confirmados_por_dia = confirmados_df[ confirmados_df['fecha_sintomas']==day.date().__str__() ] #.append(
            num_confirmados_dia.append( len(confirmados_por_dia) ) #[-1]            
            confirmados_por_dia_ingreso = confirmados_df[ confirmados_df['fecha_ingreso']==day.date().__str__() ] #.append(
            num_confirmados_por_dia_ingreso.append( len(confirmados_por_dia_ingreso) ) #[-1]
            historic_days.append( '{}/{}'.format(day.date().day,day.date().month) )
        
        fig = plt.figure(figsize=(20.0, 15.0));title="Muertes y confirmados daily historic"  
        self.plotHistoricDaily(fig, num_muertos_dia, metric, historic_days, 311 )
        self.plotHistoricDaily(fig, num_confirmados_dia, metricConfirmados, historic_days, 312 )
        self.plotHistoricDaily(fig, num_confirmados_por_dia_ingreso, metricConfirmados+" ingreso", historic_days, 313 ) #TODO: verify symptoms start
#         fig.suptitle("{}\n {}".format(title, metricRangeStr) )
#         fig.savefig(os.path.join(mobiVisuRes,title.replace(' ', '_')+".png"), bbox_inches='tight')            
        
        historic_week=[]; muertos_week=[];  muertos_week_mult=[0]; confirmados_week=[];  confirmados_week_mult=[0];
        for wk in range( int(len(num_muertos_dia)/7) ):           #TODO: verify metric in FT.com
            muertos_week.append( np.sum( num_muertos_dia[wk*7:wk*7+7] ) ); #num_muertos_dia[wk*7:wk*7+7];wk*7;wk*7+7
            confirmados_week.append( np.sum( num_confirmados_dia[wk*7:wk*7+7] ) ); 
            historic_week.append( historic_days[wk*7+3] )
            if wk>0:                
                muertos_week_mult.append(muertos_week[-1]/muertos_week[-2])
                confirmados_week_mult.append(confirmados_week[-1]/confirmados_week[-2])
        
        bad_mults=np.bitwise_or(np.isnan(muertos_week_mult), np.isinf(muertos_week_mult))
        muertos_week_mult=np.array(muertos_week_mult); muertos_week_mult[bad_mults]=0.0
        
        bad_mults=np.bitwise_or(np.isnan(confirmados_week_mult), np.isinf(confirmados_week_mult))
        confirmados_week_mult=np.array(confirmados_week_mult); confirmados_week_mult[bad_mults]=0.0
        

        fig = plt.figure(figsize=(10.0, 10.0));
        
        
        self.plotHistoricWeekly( fig, muertos_week, muertos_week_mult.tolist(), metric, historic_week, 1 )
        self.plotHistoricWeekly( fig, confirmados_week, confirmados_week_mult.tolist(), metricConfirmados, historic_week, 2 )
        
        
        
        fig.suptitle("{}\n {}".format(title, metricRangeStr) )
        fig.savefig(os.path.join(mobiVisuRes,title.replace(' ', '_')+".png"), bbox_inches='tight')

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
