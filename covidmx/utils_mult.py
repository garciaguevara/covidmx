import numpy as np
# from pylab import *

def casosPorDia(muertos_df, confirmados_df, day, num_muertos_dia, num_confirmados_dia, num_confirmados_dia_ingreso):
    muertos_por_dia = muertos_df[ muertos_df['fecha_def']==day.date().__str__() ] #.append(
    num_muertos_dia.append( len(muertos_por_dia) ) #[-1]            
    confirmados_por_dia = confirmados_df[ confirmados_df['fecha_sintomas']==day.date().__str__() ] #.append(
    num_confirmados_dia.append( len(confirmados_por_dia) ) #[-1]            
    confirmados_dia_ingreso = confirmados_df[ confirmados_df['fecha_ingreso']==day.date().__str__() ] #.append(
    num_confirmados_dia_ingreso.append( len(confirmados_dia_ingreso) ) #[-1]
    return num_muertos_dia, num_confirmados_dia, num_confirmados_dia_ingreso

def cleanMult(muertos_week_mult):
    bad_mults=np.bitwise_or(np.isnan(muertos_week_mult), np.isinf(muertos_week_mult))
    muertos_week_mult=np.array(muertos_week_mult); muertos_week_mult[bad_mults]=0.0   
#         bad_mults=np.bitwise_or(np.isnan(confirmados_week_mult), np.isinf(confirmados_week_mult))
#         confirmados_week_mult=np.array(confirmados_week_mult); confirmados_week_mult[bad_mults]=0.0    
    return muertos_week_mult.tolist()

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]    
#     N = 3
#     cumsum, moving_aves = [0], []    
#     for i, x in enumerate(mylist, 1):
#         cumsum.append(cumsum[i-1] + x)
#         if i>=N:
#             moving_ave = (cumsum[i] - cumsum[i-N])/N
#             #can do stuff with moving_ave here
#             moving_aves.append(moving_ave)    
    return ret[n - 1:] / n



def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError( "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError ("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError ( "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


