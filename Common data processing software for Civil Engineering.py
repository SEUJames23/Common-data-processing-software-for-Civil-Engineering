import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PySimpleGUI as sg
import sympy as sp
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import pandas as pd
import numpy.fft as nf
from scipy import integrate
import math
import scipy.stats as st



sg.theme('SystemDefault')

##初始页面(window0)
layout01 = [
    [sg.Button('弹塑性设计反应谱', key='-basicresponsespecture-', size=(30, 2), font=('舒体', 18), button_color=('black', 'pink'), mouseover_colors='red')],
    [sg.Text('')],
    [sg.Button('地震动弹塑性反应谱', key='-responsespecture-', size=(30, 2), font=('舒体', 18), button_color=('black', 'pink'), mouseover_colors='red')],
    [sg.Text('')],
    [sg.Button('地震波格式转换&峰值调整', key='-earthquakewaveprocessing-', size=(30, 2), font=('舒体', 18), button_color=('black', 'pink'), mouseover_colors='red')],
    [sg.Text('')],
    [sg.Button('地震波滤波处理&时程转换', key='-earthquakewavefilter-', size=(30, 2), font=('舒体', 18), button_color=('black', 'pink'), mouseover_colors='red')],
    [sg.Text('')],
    [sg.Button('IDA方法易损性分析', key='-fragility-', size=(30, 2), font=('舒体', 18), button_color=('black', 'pink'), mouseover_colors='red')]
]

layout02 = [
    [sg.Button('滞回曲线提取骨架曲线', key='-skeletoncurve-', size=(30, 2), font=('舒体', 18), button_color=('black', 'peachpuff'), mouseover_colors='red')],
    [sg.Text('')],
    [sg.Button('骨架曲线提取峰值点、极限点、等效屈服点和屈服后刚度比', key='-postyieldingstiffness-', size=(30, 2), font=('舒体', 18), button_color=('black', 'peachpuff'), mouseover_colors='red')],
    [sg.Text('')],
    [sg.Button('试验曲线降噪', key='-curvefilter-', size=(30, 2), font=('舒体', 18), button_color=('black', 'peachpuff'), mouseover_colors='red')],
]

layout0 = [
    [sg.Text("请选择需要使用的功能", text_color='red', font=('楷体', 25))],
    [sg.Text('')],
    [sg.Frame('结构动力分析相关功能', layout01, title_color='dodgerblue', font=('楷体', 18)), sg.Frame('其他数据处理功能', layout02, title_color='dodgerblue', font=('楷体', 18))]
]

window0 = sg.Window("土木工程数据常用处理软件 ----- CHEN Youjun 开发", size=(800, 800), location=(600, 80), layout=layout0, enable_close_attempted_event=True, resizable=True, element_justification='center')
win1_active = False
win2_active = False
win3_active = False
win4_active = False
win5_active = False
win6_active = False
win7_active = False
win8_active = False


while True:
    event0, values0 = window0.read()
    
    if event0 == None:
        break

## -----------window1----------window1----------window1----------window1----------window1----------window1----------window1----------
    if event0 == '-responsespecture-' and not win1_active:
        win1_active=True
        window0.Hide()
        
        ##弹塑性反应谱计算(window1)
        table_content = []
        total_filenames_location = []

        layout11 = [
            [sg.FilesBrowse("添加地震波文件", key='-openfiles-', font=20, file_types=(('ALL Files', '*.txt'), ('ALL Files', '*.AT2')), target='-filenames-'), sg.In("", key='-filenames-', font=20), sg.B('添加', key='-confirm1-', font=20), sg.B("清除选择", key="-clearchoose-", font=20), sg.B("删除", key="-delete-", font=20)],
            [sg.Table(headings=['序号', '地震名称', '时长', 'PGA'],
                    values=table_content,
                    expand_x=True,
                    expand_y=True, 
                    font=20,
                    key='-table-',
                    justification='center',
                    select_mode="extended",
                    enable_events=True)]
        ]

        layout12 = [
            [sg.Radio('加速度反应谱', group_id=1, key='加速度反应谱', default=True, font=20), sg.Radio('速度反应谱', group_id=1, key='速度反应谱', font=20), sg.Radio('位移反应谱', group_id=1, key='位移反应谱', font=20)],
            [sg.Text('阻尼比', font=20), sg.Input('0.05', key='-damping-', size=(8, 1), font=20), sg.B('确定并绘图', key='-confirm2-', font=20)],
            [sg.Canvas(key= '-toolbar-', size=(100, 1))],
            [sg.Canvas(key='-canvas-', expand_x=True, expand_y=True)],
            [sg.Text('结构周期T', font=20), sg.In('', key='-period-', size=(8, 1), font=20), sg.Text('对应平均值', font=20), sg.In('', key='-meanvalue-', size=(8, 1), font=20), sg.B('查询', key='-confirm3-', font=20)],
            [sg.FileSaveAs('创建保存文件名', target='-savelocation-', font=20, file_types=(('*.txt', '*.txt'),)), sg.In('', key='-savelocation-', visible=True, font=20), sg.B('保存', key='-save-', font=20), sg.B('清除所有', key='-clearall-', font=20), sg.B('返回至初始界面', key='-back-', button_color=('red', 'skyblue'), font=20)]
        ]

        layout10 = [
            # [sg.Menu([["帮助", ['使用说明::readme', '联系作者::contact']]], font=20)],
            [sg.Frame('地震波选择', layout11, title_color='red', font=('楷体'), expand_x=True, expand_y=True)],
            [sg.Frame('弹塑性反应谱曲线', layout12, title_color='red', font=('楷体'), expand_x=True, element_justification='center', expand_y=True)],
        ]

        layout1 = [[sg.Column(layout10, scrollable=True, vertical_scroll_only=True, size=(1200, 1100))]]

        window1 = sg.Window("土木工程数据常用处理软件 ----- 结构弹塑性反应谱计算", layout1, location=(500, 10), finalize=True, resizable=True, enable_close_attempted_event=True)

        ##win1中matplotlib初始及更新
        ##函数定义区
        def solve_sdof_eqwave_piecewise_exact(omg, zeta, ag, dt):
            omg_d = omg*np.sqrt(1.0-zeta*zeta)
            n = len(ag)
            u = np.zeros(n)
            v = np.zeros(n)
            B1 = np.exp(-zeta*omg*dt)*np.cos(omg_d*dt)
            B2 = np.exp(-zeta*omg*dt)*np.sin(omg_d*dt)
            omg_2 = 1.0/omg/omg
            omg_3 = 1.0/omg/omg/omg

            for i in range(n-1):
                u_i = u[i]
                v_i = v[i]
                p_i = -ag[i]
                alpha_i = (-ag[i+1]+ag[i])/dt

                A0 = p_i*omg_2 - 2.0*zeta*alpha_i*omg_3
                A1 = alpha_i*omg_2
                A2 = u_i - A0
                A3 = (v_i + zeta*omg*A2 - A1)/omg_d

                u[i+1] = A0 + A1*dt + A2*B1 + A3*B2
                v[i+1] = A1 + (omg_d*A3-zeta*omg*A2)*B1 - (omg_d*A2+zeta*omg*A3)*B2

            return u, v

        def response_spectra(a, dt, T, zeta):
            N = len(T)
            RSA = np.zeros(N)
            RSV = np.zeros(N)
            RSD = np.zeros(N)
            for i in range(N):
                omg = 2.0*np.pi/T[i]
                u, v = solve_sdof_eqwave_piecewise_exact(omg, zeta, ag, dt)
                a = - 2.0*zeta*omg*v - omg*omg*u
                RSA[i] = max(abs(a))
                RSV[i] = max(abs(v))
                RSD[i] = max(abs(u))
            return RSA, RSV, RSD

        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)
        fig = matplotlib.figure.Figure()
        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False
        fig.add_subplot(111)
        axes = fig.axes
        axes[0].set_xlabel('周期')
        axes[0].set_ylabel('Sa')
        figure_canvas_agg = FigureCanvasTkAgg(fig, window1['-canvas-'].TKCanvas)
        toolbar = Toolbar(figure_canvas_agg, window1['-toolbar-'].TKCanvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack()

        def updatefig(Ts, reponses, alpha, label, ylabel):
            axes=fig.axes
            axes[0].plot(Ts, reponses, alpha=alpha, label=label)
            axes[0].legend()
            axes[0].set_xlim(0, 6)
            axes[0].set_xlabel('周期')
            axes[0].set_ylabel(f'{ylabel}')
            toolbar.update()
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def scatterfig(actual_T, actual_reponse):
            axes=fig.axes
            axes[0].scatter(actual_T, actual_reponse, c='black', zorder=2)
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def clearfig():
            axes=fig.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_agg):
            figure_canvas_agg.get_tk_widget().forget()
            plt.close('all')


        while True:
            event1, values1 = window1.read()
            
            if event1 == None:
                break

            if event1 == '-confirm1-':
                try:
                    window1['-filenames-'].update('')
                    filenames = str(values1['-filenames-']).split(';')
                    current_filenames = []
                    for i in range(len(table_content)):
                        current_filenames.append(table_content[i][1])
                    for filename in filenames:
                        earthquakename = filename.split('/')[-1]

                        if earthquakename in current_filenames:
                            sg.popup('请勿重复添加地震波', title='提示', text_color='red', font=20)
                            break

                        else:
                            total_filenames_location.append(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            time = datas[0][-1]
                            PGAvalue = max(abs(datas[1]))
                            table_content.append([len(table_content)+1, earthquakename, time, PGAvalue])
                            window1['-table-'].update(table_content)
                
                except:
                    sg.popup('请先选择要打开的地震波文件', title='提示', text_color='red', font=20)

            if event1 == "-clearchoose-":
                window1["-filenames-"].update("")
            
            if event1 == "-delete-":
                index = values1["-table-"]
                num = -1
                for i in index:
                    num += 1
                    id = i - num
                    table_content.pop(id)
                    total_filenames_location.pop(id)
                for id in range(len(table_content)):
                    table_content[id][0] = id+1
                window1['-table-'].update(table_content)

            if event1 == '-confirm2-':
                clearfig()
                window1["-period-"].update("")
                window1["-meanvalue-"].update("")
                try:
                    no = 0
                    T = np.linspace(0.01, 6, 61)
                    zeta = float(values1['-damping-'])
                    RSAs, RSVs, RSDs = [], [], []
                    for filename_location in total_filenames_location:
                        no += 1
                        datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                        ag = datas[1]
                        times = datas[0]
                        dt = times[1] - times[0]
                        RSA, RSV, RSD = response_spectra(ag, dt, T, zeta)
                        RSAs.append(list(RSA))
                        RSVs.append(list(RSV))
                        RSDs.append(list(RSD))
                        if values1["加速度反应谱"]:
                            updatefig(T, RSA, 0.3, no, 'Sa')
                        if values1["速度反应谱"]:
                            updatefig(T, RSV, 0.3, no, 'Sv')
                        if values1["位移反应谱"]:
                            updatefig(T, RSD, 0.3, no, 'Sd')
                        if int(values1["加速度反应谱"])+int(values1["速度反应谱"])+int(values1["位移反应谱"])==0:
                            m=1/0 #报错用
                    
                    RSAs, RSVs, RSDs = np.array(RSAs), np.array(RSVs), np.array(RSDs)
                    mean_RSAs, mean_RSVs, mean_RSDs = RSAs.mean(axis=0), RSVs.mean(axis=0), RSDs.mean(axis=0)
                    if values1["加速度反应谱"]:
                        updatefig(T, mean_RSAs, 1, "mean", 'Sa')
                    if values1["速度反应谱"]:
                        updatefig(T, mean_RSVs, 1, "mean", 'Sv')
                    if values1["位移反应谱"]:
                        updatefig(T, mean_RSDs, 1, "mean", 'Sd')
                except ZeroDivisionError:
                    sg.popup('请检查是否选择反应谱类型', title='提示', font=20)
                except:
                    sg.popup("请检查是否添加地震波", title='提示', font=20)
            
            if event1 =="-confirm3-":
                actual_T = float(values1['-period-'])
                for id in range(len(T)-1):
                    if T[id]<=actual_T and T[id+1]>actual_T:
                        if values1["加速度反应谱"]:
                            actual_reponse = (mean_RSAs[id+1]-mean_RSAs[id])/(T[id+1]-T[id])*(actual_T-T[id])+mean_RSAs[id]
                            window1['-meanvalue-'].update('%f' %actual_reponse)
                            clearfig()
                            no = 0
                            RSAs = []
                            for filename in total_filenames_location:
                                no += 1
                                datas = np.loadtxt(filename, dtype=float, unpack=True)
                                ag = datas[1]
                                times = datas[0]
                                dt = times[1] - times[0]
                                RSA, RSV, RSD = response_spectra(ag, dt, T, zeta)
                                RSAs.append(list(RSA))
                                updatefig(T, RSA, 0.3, no, 'Sa')
                            mean_RSAs = np.array(RSAs).mean(axis=0)
                            updatefig(T, mean_RSAs, 1, "mean", 'Sa')
                            scatterfig(actual_T, actual_reponse)
                        
                        if values1["速度反应谱"]:
                            actual_reponse = (mean_RSVs[id+1]-mean_RSVs[id])/(T[id+1]-T[id])*(actual_T-T[id])+mean_RSVs[id]
                            window1['-meanvalue-'].update('%f' %actual_reponse)
                            clearfig()
                            no = 0
                            RSVs = []
                            for filename in total_filenames_location:
                                no += 1
                                datas = np.loadtxt(filename, dtype=float, unpack=True)
                                ag = datas[1]
                                times = datas[0]
                                dt = times[1] - times[0]
                                RSA, RSV, RSD = response_spectra(ag, dt, T, zeta)
                                RSVs.append(list(RSV))
                                updatefig(T, RSV, 0.3, no, 'Sv')
                            mean_RSVs = np.array(RSVs).mean(axis=0)
                            updatefig(T, mean_RSVs, 1, "mean", 'Sv')
                            scatterfig(actual_T, actual_reponse)
                        
                        if values1["位移反应谱"]:
                            actual_reponse = (mean_RSDs[id+1]-mean_RSDs[id])/(T[id+1]-T[id])*(actual_T-T[id])+mean_RSDs[id]
                            window1['-meanvalue-'].update('%f' %actual_reponse)
                            clearfig()
                            no = 0
                            RSDs = []
                            for filename in total_filenames_location:
                                no += 1
                                datas = np.loadtxt(filename, dtype=float, unpack=True)
                                ag = datas[1]
                                times = datas[0]
                                dt = times[1] - times[0]
                                RSA, RSV, RSD = response_spectra(ag, dt, T, zeta)
                                RSDs.append(list(RSD))
                                updatefig(T, RSD, 0.3, no, 'Sd')
                            mean_RSDs = np.array(RSDs).mean(axis=0)
                            updatefig(T, mean_RSDs, 1, "mean", 'Sd')
                            scatterfig(actual_T, actual_reponse)

            if event1 == "-save-":
                try:
                    savefile_location = values1['-savelocation-']
                    with open(savefile_location, 'w') as f:
                        no = 0
                        f.writelines(f"阻尼比:{zeta}\n")
                        for filename in total_filenames_location:
                            no += 1
                            earthquakename = filename.split('/')[-1]
                            f.writelines(f"地震波{no}: "+earthquakename+"\n")
                            f.writelines('周期 加速度 速度 位移\n')
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            ag = datas[1]
                            times = datas[0]
                            dt = times[1] - times[0]
                            RSA, RSV, RSD = response_spectra(ag, dt, T, zeta)
                            for t, a, v, d in zip(T, RSA, RSV, RSD):
                                f.writelines("%f %f %f %f\n" %(t, a, v, d))
                            f.writelines('\n')
                        
                        f.writelines("平均反应谱:")
                        f.writelines('周期 加速度 速度 位移\n')
                        for t, mean_a, mean_v, mean_d in zip(T, mean_RSAs, mean_RSVs, mean_RSDs):
                            f.writelines("%f %f %f %f\n" %(t, mean_a, mean_v, mean_d))
                    sg.popup_timed("保存成功！", title="提示", font=20)
                except:
                    sg.popup("操作错误!", title="提示", font=20)
            
            if event1 == '-clearall-':
                window1['-damping-'].update('0.05')
                window1['-period-'].update('')
                window1['-meanvalue-'].update('')
                window1['-savelocation-'].update('')
                window1['-filenames-'].update('')
                table_content=[]
                window1['-table-'].update(table_content)
                window1['加速度反应谱'].update(value=True)
                delete_fig_agg(figure_canvas_agg)
                total_filenames_location=[]
                sg.popup('清除成功！', title='提示', font=20)
            
            if event1 == sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
                    window1.close()
                    window0.close()
                
            if event1 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提醒', text_color='red', font=20)=="Yes":
                    window1.close()
                    win1_active=False
                    window0.UnHide()

            # if event1 == '使用说明::readme':
            #     sg.popup(image=)

            if event1 == '联系作者::contact':
                sg.popup('请发送需求到邮箱: 1765568269@qq.com', font=20)


# ## -----------window2----------window2----------window2----------window2----------window2----------window2----------window2----------
    if event0 == '-postyieldingstiffness-' and not win2_active:
        win2_active = True
        window0.Hide()
        ##等效屈服点、屈服后刚度比计算(window2)
        table_content = []
        current_filenames_location = []
        xuanze1 = ['等能量法', 'FEMA-356规定方法', '最远点法']
        xuanze2 = ['普通RC结构', '具有较高屈服后刚度比的结构']

        layout21 = [
            [sg.FilesBrowse("添加骨架曲线文件", key='-openfiles-', font=20, file_types=((('ALL Files', '*.txt'),)), target='-filenames-'), sg.In("", key='-filenames-', font=20), sg.B('添加', key='-confirm1-', font=20), sg.B("清除选择", key="-clearchoose-", font=20), sg.B("删除", key="-delete-", font=20)],
            [sg.Table(headings=['序号', '文件名称'],
                    values=table_content,
                    expand_x=True,
                    expand_y=True, 
                    font=20,
                    key='-table1-',
                    justification='center',
                    select_mode="extended",
                    enable_events=True,
                    num_rows=4)],
            [sg.Table(headings=['特征点', '位移', '荷载'], values=table_content, num_rows=4, font=20, expand_x=True, key='-table2-', justification='center')]
        ]

        layout22 = [
            [sg.Text('使用的等效方法', font=20), sg.Combo(xuanze1, key='-combo1-', font=20, default_value='等能量法', size=(18, 1)), sg.Text('结构类型', font=20), sg.Combo(xuanze2, key='-combo2-', font=20, default_value='普通RC结构', size=(26, 1)), sg.B('确定并绘图', key='-confirm2-', font=20)],
            [sg.Canvas(key='-toolbar-', size=(80, 1))], 
            [sg.Canvas(key='-canvas-', expand_x=True, expand_y=True)],
            [sg.Text("等效初始刚度", font=20), sg.In('', key='-prestiff-', disabled=True, font=20, size=(10, 1)), sg.Text("等效屈服后刚度", font=20), sg.In('', key='-poststiff-', disabled=True, font=20, size=(10, 1)), sg.Text("屈服后刚度比", font=20), sg.In('', key='-stiffratio-', disabled=True, font=20, size=(10, 1))],
            [sg.FolderBrowse('批量保存文件夹', target='-savelocation-', font=20), sg.In('', key='-savelocation-', visible=True, font=20), sg.B('保存', key='-save-', font=20), sg.B('清除所有', key='-clearall-', font=20), sg.B('返回至初始界面', key='-back-', font=20, button_color=('red', 'skyblue'))]
        ]

        layout20 = [
            [sg.Frame('结构骨架曲线选择', layout21, title_color='red', font=('楷体'), expand_x=True, expand_y=True)],
            [sg.Frame('等效双线性恢复力模型', layout22, title_color='red', font=('楷体'), expand_x=True, element_justification='center', expand_y=True)],
        ]

        layout2 = [[sg.Column(layout20, scrollable=True, vertical_scroll_only=True, size=(1200, 1100))]]
        window2 = sg.Window("土木工程数据常用处理软件 ----- 等效屈服点、屈服后刚度比计算", layout2, size=(1200, 1100), location=(300, 10), finalize=True, resizable=True, enable_close_attempted_event=True)

        def fema356method_high(filename):
            datas = np.loadtxt(filename, dtype=float, unpack=True)
            initial_disps = datas[0]
            initial_forces = datas[1]
            initial_disp0, initial_force0 = initial_disps[0], initial_forces[0]
            disps = initial_disps - np.ones_like(initial_disps)*initial_disp0
            forces = initial_forces - np.ones_like(initial_forces)*initial_force0
            k0 = forces[1]/disps[1]
            keff = 0.6*k0
            lim_f = max_f = np.max(forces)
            lim_f_d = max_f_d = float(disps[np.where(forces==max_f)][0])
            index = int(np.where(forces==max_f)[0][0])
            sum_delta_force = 0
            for id in range(index):
                id_disps = np.linspace(disps[id], disps[id+1], 100)
                id_forces = np.linspace(forces[id], forces[id+1], 100)
                delta_d = (disps[id+1] - disps[id]) / 100
                for id2 in range(len(id_disps)-1):
                    delta_force = (float(id_forces[id2]) - float(id_disps[id2])*keff) * delta_d
                    sum_delta_force += delta_force
            a = sp.symbols("a")
            tri_delta_force = (max_f_d - a) * (keff * max_f_d - max_f) / 2
            result = sp.solve([tri_delta_force + sum_delta_force], [a])
            yield_disp = result[a]
            yield_force = keff * yield_disp
            k_2 = (max_f - yield_force) / (max_f_d - yield_disp)
            ratio_k = k_2 / keff
            initial_yield_disp = yield_disp + initial_disp0
            initial_yield_force = yield_force + initial_force0
            initial_max_f_d = max_f_d + initial_disp0
            initial_max_f = max_f + initial_force0
            initial_lim_f_d = lim_f_d + initial_disp0
            initial_lim_f = lim_f + initial_force0
            return initial_disp0, initial_force0, initial_yield_disp, initial_yield_force, initial_max_f_d, initial_max_f, initial_lim_f_d, initial_lim_f, keff, k_2, ratio_k
        
        def fema356method_low(filename):
            datas = np.loadtxt(filename, dtype=float, unpack=True)
            initial_disps = datas[0]
            initial_forces = datas[1]
            initial_disp0, initial_force0 = initial_disps[0], initial_forces[0]
            disps = initial_disps - np.ones_like(initial_disps)*initial_disp0
            forces = initial_forces - np.ones_like(initial_forces)*initial_force0
            k0 = forces[1]/disps[1]
            keff = 0.6*k0
            max_f = np.max(forces)
            max_f_d = float(disps[np.where(forces==max_f)][0])
            index = int(np.where(forces==max_f)[0][0])
            sum_delta_force = 0
            for id in range(index):
                id_disps = np.linspace(disps[id], disps[id+1], 100)
                id_forces = np.linspace(forces[id], forces[id+1], 100)
                delta_d = (disps[id+1] - disps[id]) / 100
                for id2 in range(len(id_disps)-1):
                    delta_force = (float(id_forces[id2]) - float(id_disps[id2])*keff) * delta_d
                    sum_delta_force += delta_force
            a = sp.symbols("a")
            tri_delta_force = (max_f_d - a) * (keff * max_f_d - max_f) / 2
            result = sp.solve([tri_delta_force + sum_delta_force], [a])
            yield_disp = result[a]
            yield_force = keff * yield_disp
            lim_f = 0.85 * max_f
            for no in range(index, len(disps)):
                if forces[no]>=lim_f and forces[no+1]<lim_f:
                    lim_f_d = (disps[no+1] - disps[no]) / (forces[no+1] - forces[no]) * (lim_f - forces[no]) + disps[no]
            k_2 = (lim_f - yield_force) / (lim_f_d - yield_disp)
            ratio_k = k_2 / keff
            initial_yield_disp = yield_disp + initial_disp0
            initial_yield_force = yield_force + initial_force0
            initial_max_f_d = max_f_d + initial_disp0
            initial_max_f = max_f + initial_force0
            initial_lim_f_d = lim_f_d + initial_disp0
            initial_lim_f = lim_f + initial_force0
            return initial_disp0, initial_force0, initial_yield_disp, initial_yield_force, initial_max_f_d, initial_max_f, initial_lim_f_d, initial_lim_f, keff, k_2, ratio_k
        
        def equalenergymethod_high(filename):
            datas = np.loadtxt(filename, dtype=float, unpack=True)
            initial_disps = datas[0]
            initial_forces = datas[1]
            initial_disp0, initial_force0 = initial_disps[0], initial_forces[0]
            disps = initial_disps - np.ones_like(initial_disps)*initial_disp0
            forces = initial_forces - np.ones_like(initial_forces)*initial_force0
            lim_f = max_f = np.max(forces)
            lim_f_d = max_f_d = float(disps[np.where(forces==max_f)][0])
            index = int(np.where(forces==max_f)[0][0])
            new_disps = disps[:index+1]
            new_forces = forces[:index+1]
            sum_delta_force = 0
            for id in range(len(new_disps)-1):
                id_disps = np.linspace(new_disps[id], new_disps[id+1], 100)
                id_forces = np.linspace(new_forces[id], new_forces[id+1], 100)
                delta_d = (disps[id+1] - disps[id]) / 100
                for id2 in range(len(id_disps)-1):
                    delta_force = (float(id_forces[id2]) + float(id_forces[id2+1])) * delta_d / 2
                    sum_delta_force += delta_force
            a = sp.symbols("a")
            area = (2 * max_f_d - a) * max_f / 2
            result = sp.solve([area - sum_delta_force], [a])
            keff = max_f/result[a]
            for id in range(len(disps)):
                if forces[id]>=float(keff*disps[id]) and forces[id+1]<float(keff*disps[id+1]):
                    b = sp.symbols("b")
                    km = (forces[id+1]-forces[id])/(disps[id+1]-disps[id])
                    mm = forces[id]-km*disps[id]
                    met1 = km * b + mm
                    met2 = keff * b
                    result = sp.solve([met1-met2], [b])
                    yield_disp = result[b]
            yield_force = yield_disp*keff
            k_2 = (max_f - yield_force) / (max_f_d - yield_disp)
            ratio_k = k_2 / keff
            initial_yield_disp = yield_disp + initial_disp0
            initial_yield_force = yield_force + initial_force0
            initial_max_f_d = max_f_d + initial_disp0
            initial_max_f = max_f + initial_force0
            initial_lim_f_d = lim_f_d + initial_disp0
            initial_lim_f = lim_f + initial_force0
            return initial_disp0, initial_force0, initial_yield_disp, initial_yield_force, initial_max_f_d, initial_max_f, initial_lim_f_d, initial_lim_f, keff, k_2, ratio_k
        
        
        def equalenergymethod_low(filename):
            datas = np.loadtxt(filename, dtype=float, unpack=True)
            initial_disps = datas[0]
            initial_forces = datas[1]
            initial_disp0, initial_force0 = initial_disps[0], initial_forces[0]
            disps = initial_disps - np.ones_like(initial_disps)*initial_disp0
            forces = initial_forces - np.ones_like(initial_forces)*initial_force0
            
            max_f = np.max(forces)
            max_f_d = float(disps[np.where(forces==max_f)][0])
            index = int(np.where(forces==max_f)[0][0])
            new_disps = disps[:index+1]
            new_forces = forces[:index+1]
            sum_delta_force = 0
            for id in range(len(new_disps)-1):
                id_disps = np.linspace(new_disps[id], new_disps[id+1], 100)
                id_forces = np.linspace(new_forces[id], new_forces[id+1], 100)
                delta_d = (disps[id+1] - disps[id]) / 100
                for id2 in range(len(id_disps)-1):
                    delta_force = (float(id_forces[id2]) + float(id_forces[id2+1])) * delta_d / 2
                    sum_delta_force += delta_force
            
            a = sp.symbols("a")
            area = (2 * max_f_d - a) * max_f / 2
            result = sp.solve([area - sum_delta_force], [a])
            keff = max_f/result[a]
            for id in range(len(disps)):
                if forces[id]>=float(keff*disps[id]) and forces[id+1]<float(keff*disps[id+1]):
                    b = sp.symbols("b")
                    km = (forces[id+1]-forces[id])/(disps[id+1]-disps[id])
                    mm = forces[id]-km*disps[id]
                    met1 = km * b + mm
                    met2 = keff * b
                    result = sp.solve([met1-met2], [b])
                    yield_disp = result[b]
            yield_force = yield_disp*keff
            lim_f = 0.85 * max_f
            for no in range(index, len(disps)):
                if forces[no]>=lim_f and forces[no+1]<lim_f:
                    lim_f_d = (disps[no+1] - disps[no]) / (forces[no+1] - forces[no]) * (lim_f - forces[no]) + disps[no]
            k_2 = (lim_f - max_f) / (lim_f_d - max_f_d)
            ratio_k = k_2 / keff
            initial_yield_disp = yield_disp + initial_disp0
            initial_yield_force = yield_force + initial_force0
            initial_max_f_d = max_f_d + initial_disp0
            initial_max_f = max_f + initial_force0
            initial_lim_f_d = lim_f_d + initial_disp0
            initial_lim_f = lim_f + initial_force0
            return initial_disp0, initial_force0, initial_yield_disp, initial_yield_force, initial_max_f_d, initial_max_f, initial_lim_f_d, initial_lim_f, keff, k_2, ratio_k
        
        
        def farthestpointmethod_high(filename):
            datas = np.loadtxt(filename, dtype=float, unpack=True)
            initial_disps = datas[0]
            initial_forces = datas[1]
            initial_disp0, initial_force0 = initial_disps[0], initial_forces[0]
            disps = initial_disps - np.ones_like(initial_disps)*initial_disp0
            forces = initial_forces - np.ones_like(initial_forces)*initial_force0
            lim_f = max_f = np.max(forces)
            lim_f_d = max_f_d = float(disps[np.where(forces==max_f)][0])
            index = int(np.where(forces==max_f)[0][0])
            k1 = max_f/max_f_d
            new_disps = disps[:index+1]
            new_forces = forces[:index+1]
            divide_disps, divide_forces = [], []
            for id in range(len(new_disps)-1):
                id_disps = np.linspace(new_disps[id], new_disps[id+1], 100)
                id_forces = np.linspace(new_forces[id], new_forces[id+1], 100)
                for disp0, force0 in zip(id_disps, id_forces):
                    divide_disps.append(float(disp0))
                    divide_forces.append(float(force0))
            distance_fenzi = []
            for di_disp, di_force in zip(divide_disps, divide_forces):
                dist = abs(k1*di_disp-di_force+max_f-k1*max_f_d)
                distance_fenzi.append(dist)
            distance_fenzi = np.array(distance_fenzi).reshape(-1, 1)
            distance_max = np.max(distance_fenzi)
            farthestid = 0
            for id_max in range(len(distance_fenzi)):
                if distance_fenzi[id_max]==distance_max:
                    farthestid = id_max
            yield_disp = divide_disps[farthestid]
            yield_force = divide_forces[farthestid]
            keff = yield_force/yield_disp
            k_2 = (max_f-yield_force)/(max_f_d-yield_disp)
            ratio_k = k_2 / keff
            initial_yield_disp = yield_disp + initial_disp0
            initial_yield_force = yield_force + initial_force0
            initial_max_f_d = max_f_d + initial_disp0
            initial_max_f = max_f + initial_force0
            initial_lim_f_d = lim_f_d + initial_disp0
            initial_lim_f = lim_f + initial_force0
            return initial_disp0, initial_force0, initial_yield_disp, initial_yield_force, initial_max_f_d, initial_max_f, initial_lim_f_d, initial_lim_f, keff, k_2, ratio_k
        
        
        def farthestpointmethod_low(filename):
            datas = np.loadtxt(filename, dtype=float, unpack=True)
            initial_disps = datas[0]
            initial_forces = datas[1]
            initial_disp0, initial_force0 = initial_disps[0], initial_forces[0]
            disps = initial_disps - np.ones_like(initial_disps)*initial_disp0
            forces = initial_forces - np.ones_like(initial_forces)*initial_force0
            max_f = np.max(forces)
            max_f_d = float(disps[np.where(forces==max_f)][0])
            index = int(np.where(forces==max_f)[0][0])
            lim_f = 0.85 * max_f
            for no in range(index, len(disps)):
                if forces[no]>=lim_f and forces[no+1]<lim_f:
                    lim_f_d = (disps[no+1] - disps[no]) / (forces[no+1] - forces[no]) * (lim_f - forces[no]) + disps[no]
            k1 = lim_f/lim_f_d
            new_disps = disps[:index+1]
            new_forces = forces[:index+1]
            divide_disps, divide_forces = [], []
            for id in range(len(new_disps)-1):
                id_disps = np.linspace(new_disps[id], new_disps[id+1], 100)
                id_forces = np.linspace(new_forces[id], new_forces[id+1], 100)
                for disp0, force0 in zip(id_disps, id_forces):
                    divide_disps.append(float(disp0))
                    divide_forces.append(float(force0))
            distance_fenzi = []
            for di_disp, di_force in zip(divide_disps, divide_forces):
                dist = abs(k1*di_disp-di_force+lim_f-k1*lim_f_d)
                distance_fenzi.append(dist)
            distance_fenzi = np.array(distance_fenzi).reshape(-1, 1)
            distance_max = np.max(distance_fenzi)
            farthestid = 0
            for id_max in range(len(distance_fenzi)):
                if distance_fenzi[id_max]==distance_max:
                    farthestid = id_max
            yield_disp = divide_disps[farthestid]
            yield_force = divide_forces[farthestid]
            keff = yield_force/yield_disp
            k_2 = (lim_f-yield_force)/(lim_f_d-yield_disp)
            ratio_k = k_2 / keff
            initial_yield_disp = yield_disp + initial_disp0
            initial_yield_force = yield_force + initial_force0
            initial_max_f_d = max_f_d + initial_disp0
            initial_max_f = max_f + initial_force0
            initial_lim_f_d = lim_f_d + initial_disp0
            initial_lim_f = lim_f + initial_force0
            return initial_disp0, initial_force0, initial_yield_disp, initial_yield_force, initial_max_f_d, initial_max_f, initial_lim_f_d, initial_lim_f, keff, k_2, ratio_k
        
        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)

        fig = matplotlib.figure.Figure(figsize=(8, 5.5))
        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False
        fig.add_subplot(111)
        axes = fig.axes
        axes[0].set_xlabel('位移')
        axes[0].set_ylabel('荷载')
        figure_canvas_agg = FigureCanvasTkAgg(fig, window2['-canvas-'].TKCanvas)
        toolbar = Toolbar(figure_canvas_agg, window2['-toolbar-'].TKCanvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack()

        def updatefig(disps, forces, alpha, label):
            axes=fig.axes
            axes[0].set_xlabel('位移')
            axes[0].set_ylabel('荷载')
            axes[0].plot(disps, forces, alpha=alpha, label=label)
            axes[0].legend()
            axes[0].set_xlim(0, )
            toolbar.update()
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def scatterfig(point_disp, point_force):
            axes=fig.axes
            axes[0].scatter(point_disp, point_force, c='black', zorder=2)
            axes[0].set_ylim(0, )
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def clearfig():
            axes=fig.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_agg):
            figure_canvas_agg.get_tk_widget().forget()
            plt.close('all')

        while True:
            event2, values2 = window2.read()

            if event2 == None:
                break

            if event2 == '-confirm1-':
                try:
                    window2['-filenames-'].update('')
                    filenames = str(values2['-filenames-']).split(';')
                    current_filenames = []
                    for i in range(len(table_content)):
                        current_filenames.append(table_content[i][1])
                    for filename in filenames:
                        dataname = filename.split('/')[-1]
                        if dataname in current_filenames:
                            sg.popup('请勿重复添加文件', title='提示', text_color='red', font=20)
                            break
                        else:
                            current_filenames_location.append(filename)
                            table_content.append([len(table_content)+1, dataname])
                            window2['-table1-'].update(table_content)
                except:
                    sg.popup('请先选择要打开的文件', title='提示', text_color='red', font=20)
            
            if event2 == "-clearchoose-":
                window2["-filenames-"].update("")
            
            if event2 == "-delete-":
                index = values2["-table1-"]
                num = -1
                for i in index:
                    num += 1
                    id = i - num
                    table_content.pop(id)
                    current_filenames_location.pop(id)
                for id in range(len(table_content)):
                    table_content[id][0] = id+1
                window2['-table1-'].update(table_content)

            if event2 == '-confirm2-':
                clearfig()
                table_content2=[]
                window2['-table2-'].update(table_content2)
                try:
                    if values2['-combo1-']=='等能量法' and values2['-combo2-']=='普通RC结构':
                        index = values2['-table1-']
                        if len(index)>1:
                            sg.popup("图中仅限一组数据的骨架曲线,请在上表中仅选择一组", title='提示', text_color='red', font=20)
                        elif len(index)==0:
                            sg.popup("请在上表中选择一组需要绘图的数据", title='提示', text_color='red', font=20)
                        else:
                            filename = current_filenames_location[index[0]]
                            initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = equalenergymethod_low(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            disps = datas[0]
                            forces = datas[1]
                            updatefig(disps, forces, 1, label='骨架曲线')
                            updatefig([initial_disp, yield_disp], [initial_force, yield_force], 1, label='等效初始刚度')
                            updatefig([max_f_d, lim_f_d], [max_f, lim_f], 1, label='等效屈服后刚度')
                            scatterfig(initial_disp, initial_force)
                            scatterfig(yield_disp, yield_force)
                            scatterfig(max_f_d, max_f)
                            scatterfig(lim_f_d, lim_f)
                            table_content2 = [['初始点', str(round(initial_disp, 2)), str(round(initial_force, 2))], ['等效屈服点', str(round(yield_disp, 2)), str(round(yield_force, 2))], ['峰值点', str(round(max_f_d, 2)), str(round(max_f, 2))], ['极限点', str(round(lim_f_d, 2)), str(round(lim_f, 2))]]
                            window2['-table2-'].update(table_content2)
                            window2['-prestiff-'].update(f'{round(keff, 3)}')
                            window2['-poststiff-'].update(f'{round(k_2, 3)}')
                            window2['-stiffratio-'].update(f'{round(ratio_k, 4)}')
                    
                    if values2['-combo1-']=='等能量法' and values2['-combo2-']=='具有较高屈服后刚度比的结构':
                        index = values2['-table1-']
                        if len(index)>1:
                            sg.popup("图中仅限一组数据的骨架曲线,请在上表中仅选择一组", title='提示', text_color='red', font=20)
                        elif len(index)==0:
                            sg.popup("请在上表中选择一组需要绘图的数据", title='提示', text_color='red', font=20)
                        else:
                            filename = current_filenames_location[index[0]]
                            initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = equalenergymethod_high(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            disps = datas[0]
                            forces = datas[1]
                            updatefig(disps, forces, 1, label='骨架曲线')
                            updatefig([initial_disp, yield_disp, lim_f_d], [initial_force, yield_force, lim_f], 1, label='等效双线性')
                            scatterfig(initial_disp, initial_force)
                            scatterfig(yield_disp, yield_force)
                            scatterfig(lim_f_d, lim_f)
                            table_content2 = [['初始点', str(round(initial_disp, 2)), str(round(initial_force, 2))], ['等效屈服点', str(round(yield_disp, 2)), str(round(yield_force, 2))], ['峰值点', str(round(max_f_d, 2)), str(round(max_f, 2))], ['极限点', str(round(lim_f_d, 2)), str(round(lim_f, 2))]]
                            window2['-table2-'].update(table_content2)
                            window2['-prestiff-'].update(f'{round(keff, 3)}')
                            window2['-poststiff-'].update(f'{round(k_2, 3)}')
                            window2['-stiffratio-'].update(f'{round(ratio_k, 4)}')

                    if values2['-combo1-']=='FEMA-356规定方法' and values2['-combo2-']=='普通RC结构':
                        index = values2['-table1-']
                        if len(index)>1:
                            sg.popup("图中仅限一组数据的骨架曲线,请在上表中仅选择一组", title='提示', text_color='red', font=20)
                        elif len(index)==0:
                            sg.popup("请在上表中选择一组需要绘图的数据", title='提示', text_color='red', font=20)
                        else:
                            filename = current_filenames_location[index[0]]
                            initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = fema356method_low(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            disps = datas[0]
                            forces = datas[1]
                            updatefig(disps, forces, 1, label='骨架曲线')
                            updatefig([initial_disp, yield_disp, lim_f_d], [initial_force, yield_force, lim_f], 1, label='等效双线性')
                            scatterfig(initial_disp, initial_force)
                            scatterfig(yield_disp, yield_force)
                            scatterfig(lim_f_d, lim_f)
                            scatterfig(max_f_d, max_f)
                            table_content2 = [['初始点', str(round(initial_disp, 2)), str(round(initial_force, 2))], ['等效屈服点', str(round(yield_disp, 2)), str(round(yield_force, 2))], ['峰值点', str(round(max_f_d, 2)), str(round(max_f, 2))], ['极限点', str(round(lim_f_d, 2)), str(round(lim_f, 2))]]
                            window2['-table2-'].update(table_content2)
                            window2['-prestiff-'].update(f'{round(keff, 3)}')
                            window2['-poststiff-'].update(f'{round(k_2, 3)}')
                            window2['-stiffratio-'].update(f'{round(ratio_k, 4)}')

                    if values2['-combo1-']=='FEMA-356规定方法' and values2['-combo2-']=='具有较高屈服后刚度比的结构':
                        index = values2['-table1-']
                        if len(index)>1:
                            sg.popup("图中仅限一组数据的骨架曲线,请在上表中仅选择一组", title='提示', text_color='red', font=20)
                        elif len(index)==0:
                            sg.popup("请在上表中选择一组需要绘图的数据", title='提示', text_color='red', font=20)
                        else:
                            filename = current_filenames_location[index[0]]
                            initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = fema356method_high(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            disps = datas[0]
                            forces = datas[1]
                            updatefig(disps, forces, 1, label='骨架曲线')
                            updatefig([initial_disp, yield_disp, lim_f_d], [initial_force, yield_force, lim_f], 1, label='等效双线性')
                            scatterfig(initial_disp, initial_force)
                            scatterfig(yield_disp, yield_force)
                            scatterfig(lim_f_d, lim_f)
                            table_content2 = [['初始点', str(round(initial_disp, 2)), str(round(initial_force, 2))], ['等效屈服点', str(round(yield_disp, 2)), str(round(yield_force, 2))], ['峰值点', str(round(max_f_d, 2)), str(round(max_f, 2))], ['极限点', str(round(lim_f_d, 2)), str(round(lim_f, 2))]]
                            window2['-table2-'].update(table_content2)
                            window2['-prestiff-'].update(f'{round(keff, 3)}')
                            window2['-poststiff-'].update(f'{round(k_2, 3)}')
                            window2['-stiffratio-'].update(f'{round(ratio_k, 4)}')

                    if values2['-combo1-']=='最远点法' and values2['-combo2-']=='普通RC结构':
                        index = values2['-table1-']
                        if len(index)>1:
                            sg.popup("图中仅限一组数据的骨架曲线,请在上表中仅选择一组", title='提示', text_color='red', font=20)
                        elif len(index)==0:
                            sg.popup("请在上表中选择一组需要绘图的数据", title='提示', text_color='red', font=20)
                        else:
                            filename = current_filenames_location[index[0]]
                            initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = farthestpointmethod_low(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            disps = datas[0]
                            forces = datas[1]
                            updatefig(disps, forces, 1, label='骨架曲线')
                            updatefig([initial_disp, yield_disp, lim_f_d], [initial_force, yield_force, lim_f], 1, label='等效双线性')
                            scatterfig(initial_disp, initial_force)
                            scatterfig(yield_disp, yield_force)
                            scatterfig(max_f_d, max_f)
                            scatterfig(lim_f_d, lim_f)
                            table_content2 = [['初始点', str(round(initial_disp, 2)), str(round(initial_force, 2))], ['等效屈服点', str(round(yield_disp, 2)), str(round(yield_force, 2))], ['峰值点', str(round(max_f_d, 2)), str(round(max_f, 2))], ['极限点', str(round(lim_f_d, 2)), str(round(lim_f, 2))]]
                            window2['-table2-'].update(table_content2)
                            window2['-prestiff-'].update(f'{round(keff, 3)}')
                            window2['-poststiff-'].update(f'{round(k_2, 3)}')
                            window2['-stiffratio-'].update(f'{round(ratio_k, 4)}')

                    if values2['-combo1-']=='最远点法' and values2['-combo2-']=='具有较高屈服后刚度比的结构':
                        index = values2['-table1-']
                        if len(index)>1:
                            sg.popup("图中仅限一组数据的骨架曲线,请在上表中仅选择一组", title='提示', text_color='red', font=20)
                        elif len(index)==0:
                            sg.popup("请在上表中选择一组需要绘图的数据", title='提示', text_color='red', font=20)
                        else:
                            filename = current_filenames_location[index[0]]
                            initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = farthestpointmethod_high(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            disps = datas[0]
                            forces = datas[1]
                            updatefig(disps, forces, 1, label='骨架曲线')
                            updatefig([initial_disp, yield_disp, lim_f_d], [initial_force, yield_force, lim_f], 1, label='等效双线性')
                            scatterfig(initial_disp, initial_force)
                            scatterfig(yield_disp, yield_force)
                            scatterfig(max_f_d, max_f)
                            table_content2 = [['初始点', str(round(initial_disp, 2)), str(round(initial_force, 2))], ['等效屈服点', str(round(yield_disp, 2)), str(round(yield_force, 2))], ['峰值点', str(round(max_f_d, 2)), str(round(max_f, 2))], ['极限点', str(round(lim_f_d, 2)), str(round(lim_f, 2))]]
                            window2['-table2-'].update(table_content2)
                            window2['-prestiff-'].update(f'{round(keff, 3)}')
                            window2['-poststiff-'].update(f'{round(k_2, 3)}')
                            window2['-stiffratio-'].update(f'{round(ratio_k, 4)}')
                except:
                    sg.popup('请检查结构类型或者数据长度！', title='提示', text_color='red', font=20)   

            if event2=="-save-":
                try:
                    folder_location=values2["-savelocation-"]
                    index = values2['-table1-']

                    if len(index)!=0:
                        choose_filenames_location = []
                        choose_filenames = []
                        for ind in index:
                            choose_filenames_location.append(current_filenames_location[ind])
                            choose_filenames.append(table_content[ind][1])
                        
                        no = 0
                        for filename_location, filename in zip(choose_filenames_location, choose_filenames):
                            no += 1
                            savefile_location = folder_location+'/特征点提取_'+filename
                            with open(savefile_location, "w") as f:
                                if values2['-combo1-']=='等能量法' and values2['-combo2-']=='普通RC结构':
                                    initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = equalenergymethod_low(filename_location)
                                    f.writelines(f'结构{no}:'+table_content[index[no-1]][1]+"\n")
                                    f.writelines(f"使用的等效方法：{values2['-combo1-']}   结构类型：{values2['-combo2-']}\n")
                                    f.writelines(f'初始位移:{initial_disp} 初始荷载:{initial_force}\n')
                                    f.writelines(f'等效屈服位移:{yield_disp} 等效屈服荷载:{yield_force}\n')
                                    f.writelines(f'等效屈服位移:{yield_disp} 等效屈服荷载:{yield_force}\n')
                                    f.writelines(f'峰值荷载对应的位移:{max_f_d} 峰值荷载:{max_f}\n')
                                    f.writelines(f'极限荷载对应的位移:{lim_f_d} 极限荷载:{lim_f}\n')
                                    f.writelines(f'等效初始刚度:{keff} 等效屈服后刚度:{k_2}\n')
                                    f.writelines(f'屈服后刚度比:{ratio_k}\n')
                                
                                if values2['-combo1-']=='等能量法' and values2['-combo2-']=='具有较高屈服后刚度比的结构':
                                    initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = equalenergymethod_high(filename_location)
                                    f.writelines(f'结构{no}:'+table_content[index[no-1]][1]+"\n")
                                    f.writelines(f"使用的等效方法：{values2['-combo1-']}   结构类型：{values2['-combo2-']}\n")
                                    f.writelines(f'初始位移:{initial_disp} 初始荷载:{initial_force}\n')
                                    f.writelines(f'等效屈服位移:{yield_disp} 等效屈服荷载:{yield_force}\n')
                                    f.writelines(f'峰值荷载对应的位移:{max_f_d} 峰值荷载:{max_f}\n')
                                    f.writelines(f'极限荷载对应的位移:{lim_f_d} 极限荷载:{lim_f}\n')
                                    f.writelines(f'等效初始刚度:{keff} 等效屈服后刚度:{k_2}\n')
                                    f.writelines(f'屈服后刚度比:{ratio_k}\n')

                                if values2['-combo1-']=='FEMA-356规定方法' and values2['-combo2-']=='普通RC结构':
                                    initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = fema356method_low(filename_location)
                                    f.writelines(f'结构{no}:'+table_content[index[no-1]][1]+"\n")
                                    f.writelines(f"使用的等效方法：{values2['-combo1-']}   结构类型：{values2['-combo2-']}\n")
                                    f.writelines(f'初始位移:{initial_disp} 初始荷载:{initial_force}\n')
                                    f.writelines(f'等效屈服位移:{yield_disp} 等效屈服荷载:{yield_force}\n')
                                    f.writelines(f'峰值荷载对应的位移:{max_f_d} 峰值荷载:{max_f}\n')
                                    f.writelines(f'极限荷载对应的位移:{lim_f_d} 极限荷载:{lim_f}\n')
                                    f.writelines(f'等效初始刚度:{keff} 等效屈服后刚度:{k_2}\n')
                                    f.writelines(f'屈服后刚度比:{ratio_k}\n')

                                if values2['-combo1-']=='FEMA-356规定方法' and values2['-combo2-']=='具有较高屈服后刚度比的结构':
                                    initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = fema356method_high(filename_location)
                                    f.writelines(f'结构{no}:'+table_content[index[no-1]][1]+"\n")
                                    f.writelines(f"使用的等效方法：{values2['-combo1-']}   结构类型：{values2['-combo2-']}\n")
                                    f.writelines(f'初始位移:{initial_disp} 初始荷载:{initial_force}\n')
                                    f.writelines(f'等效屈服位移:{yield_disp} 等效屈服荷载:{yield_force}\n')
                                    f.writelines(f'峰值荷载对应的位移:{max_f_d} 峰值荷载:{max_f}\n')
                                    f.writelines(f'极限荷载对应的位移:{lim_f_d} 极限荷载:{lim_f}\n')
                                    f.writelines(f'等效初始刚度:{keff} 等效屈服后刚度:{k_2}\n')
                                    f.writelines(f'屈服后刚度比:{ratio_k}\n')

                                if values2['-combo1-']=='最远点法' and values2['-combo2-']=='普通RC结构':
                                    initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = farthestpointmethod_low(filename_location)
                                    f.writelines(f'结构{no}:'+table_content[index[no-1]][1]+"\n")
                                    f.writelines(f"使用的等效方法：{values2['-combo1-']}   结构类型：{values2['-combo2-']}\n")
                                    f.writelines(f'初始位移:{initial_disp} 初始荷载:{initial_force}\n')
                                    f.writelines(f'等效屈服位移:{yield_disp} 等效屈服荷载:{yield_force}\n')
                                    f.writelines(f'峰值荷载对应的位移:{max_f_d} 峰值荷载:{max_f}\n')
                                    f.writelines(f'极限荷载对应的位移:{lim_f_d} 极限荷载:{lim_f}\n')
                                    f.writelines(f'等效初始刚度:{keff} 等效屈服后刚度:{k_2}\n')
                                    f.writelines(f'屈服后刚度比:{ratio_k}\n')

                                if values2['-combo1-']=='最远点法' and values2['-combo2-']=='具有较高屈服后刚度比的结构':
                                    initial_disp, initial_force, yield_disp, yield_force, max_f_d, max_f, lim_f_d, lim_f, keff, k_2, ratio_k = farthestpointmethod_high(filename_location)
                                    f.writelines(f'结构{no}:'+table_content[index[no-1]][1]+"\n")
                                    f.writelines(f"使用的等效方法：{values2['-combo1-']}   结构类型：{values2['-combo2-']}\n")
                                    f.writelines(f'初始位移:{initial_disp} 初始荷载:{initial_force}\n')
                                    f.writelines(f'等效屈服位移:{yield_disp} 等效屈服荷载:{yield_force}\n')
                                    f.writelines(f'峰值荷载对应的位移:{max_f_d} 峰值荷载:{max_f}\n')
                                    f.writelines(f'极限荷载对应的位移:{lim_f_d} 极限荷载:{lim_f}\n')
                                    f.writelines(f'等效初始刚度:{keff} 等效屈服后刚度:{k_2}\n')
                                    f.writelines(f'屈服后刚度比:{ratio_k}\n')
                        sg.popup('保存成功！', title='提示', font=20)
                    else:
                        sg.popup('请先在上表中选择需要导出的文件名，注意区分结构类型', title='提示', text_color='red', font=20)
                except:
                    sg.popup('请检查结构类型、数据长度或者是否添加保存文件夹！', title='提示', text_color='red', font=20)       

            if event2=='-clearall-':
                window2['-filenames-'].update('')
                table_content = table_content2 = []
                window2['-table1-'].update(table_content)
                window2['-table2-'].update(table_content2)
                window2['-savelocation-'].update('')
                window2['-combo1-'].update('等能量法')
                window2['-combo2-'].update('普通RC结构')
                window2['-prestiff-'].update('')
                window2['-poststiff-'].update('')
                window2['-stiffratio-'].update('')
                delete_fig_agg(figure_canvas_agg)
                current_filenames_location=[]
                choose_filenames_location=[]
                sg.popup('清除成功！', font=20)
            
            if event2==sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
                    window2.close()
                    window0.close()

            if event2 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提醒', text_color='red', font=20)=="Yes":
                    window2.close()
                    win2_active=False
                    window0.UnHide()



## -----------window3----------window3----------window3----------window3----------window3----------window3----------window3----------
    if event0 == '-basicresponsespecture-' and not win3_active:
        win3_active=True
        window0.Hide()
        lieduxz = ['6度', '7度', '7度半', '8度', '8度半', '9度']
        shefangxz = ['多遇地震', '设防地震', '罕遇地震', '极罕遇地震']
        dizhenfenzuxz = ['第一组', '第二组', '第三组']
        changdixz = ['I0', 'I1', 'Ⅱ', 'Ⅲ', 'Ⅳ']

        layout31 = [
            [sg.Text('地震烈度', font=20), sg.Combo(lieduxz, key='-combo1-', default_value='6度', size=(15, 1), font=20), sg.Text('    '), sg.Text('设防等级', font=20), sg.Combo(shefangxz, key='-combo2-', font=20, default_value='设防地震', size=(15, 1))],
            [sg.Text('地震分组', font=20), sg.Combo(dizhenfenzuxz, key='-combo3-', font=20, default_value='第一组', size=(15, 1)), sg.Text('    '), sg.Text('场地类别', font=20), sg.Combo(changdixz, key='-combo4-', font=20, default_value='Ⅱ', size=(15, 1))],
            [sg.Text('阻尼比', font=20), sg.Input('0.05', key='-damping-', size=(15, 1), font=20)],
            [sg.B('确定并绘图', key='-confirm1-', font=20), sg.Text('    '), sg.B('清除选择', key='-delete-', font=20)],
            [sg.FileSaveAs('创建保存文件名', font=20, target='-savelocation-', file_types=(('*.txt', '*.txt'),)), sg.In('', key='-savelocation-', visible=True, font=20), sg.B('保存', key='-save-', font=20)]
        ]

        layout32 = [
            [sg.Canvas(key='-toolbar-', size=(80, 1))],
            [sg.Canvas(key='-canvas-', expand_x=True, expand_y=True)],
            [sg.Text('结构周期', font=20), sg.In('', key='-period-', size=(15, 1)), sg.Text('  '), sg.Text('对应谱加速度值', font=20), sg.In('', key='-sa-', size=(15, 1)), sg.Text('    '), sg.B('查询', key='-confirm2-', font=20), sg.Text('    '), sg.B('返回至初始界面', key='-back-', font=20, button_color=('red', 'skyblue'))]
        ]

        layout30 = [
            [sg.Frame('地震参数选择', layout31, title_color='red', font=('楷体', 15), element_justification='center', expand_x=True, expand_y=True)], 
            [sg.Frame('绘图', layout32, font=('楷体', 15), title_color='red', expand_x=True, expand_y=True, element_justification='center')]
        ]

        layout3 = [[sg.Column(layout30, scrollable=True, vertical_scroll_only=True, size=(900, 900))]]
        window3 = sg.Window("土木工程数据常用处理软件 ----- 抗震规范中弹塑性设计反应谱计算", layout3, location=(600, 120), finalize=True, resizable=True, enable_close_attempted_event=True)

        def figdata(afa_max, tg, gama, eta1, eta2):
            times, afas = np.array([]), np.array([])
            #上升段
            x1 = np.linspace(0, 0.1, 10)
            y1 = np.linspace(0.45*afa_max, eta2*afa_max, 10)
            times = np.concatenate((times, x1), axis=0)
            afas = np.concatenate((afas, y1), axis=0)
            #平台段
            x2 = np.linspace(0.1, tg, 10)
            y2 = np.ones_like(x2)*afa_max*eta2
            times = np.concatenate((times, x2), axis=0)
            afas = np.concatenate((afas, y2), axis=0)
            #曲线下降段
            x3 = np.linspace(tg, 5*tg, 100)
            y3 = pow(tg/x3, gama)*eta2*afa_max
            times = np.concatenate((times, x3), axis=0)
            afas = np.concatenate((afas, y3), axis=0)
            #直线下降段
            x4 = np.linspace(5*tg, 6, 100)
            y4 = (eta2*pow(0.2, gama)-eta1*(x4-5*tg))*afa_max
            times = np.concatenate((times, x4), axis=0)
            afas = np.concatenate((afas, y4), axis=0)
            return times, afas

        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)
        fig = matplotlib.figure.Figure()
        fig = matplotlib.figure.Figure()
        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False
        fig.add_subplot(111)
        axes = fig.axes
        axes[0].set_xlabel('周期')
        axes[0].set_ylabel('Sa')
        figure_canvas_agg = FigureCanvasTkAgg(fig, window3['-canvas-'].TKCanvas)
        toolbar = Toolbar(figure_canvas_agg, window3['-toolbar-'].TKCanvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack()

        def updatefig(xs, ys, alpha, label):
            axes=fig.axes
            axes[0].set_xlabel('周期')
            axes[0].set_ylabel('Sa')
            axes[0].plot(xs, ys, alpha=alpha, label=label)
            axes[0].set_xlim(0, 6)
            axes[0].set_ylim(0, )
            axes[0].legend()
            toolbar.update()
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def scatterfig(actual_T, actual_afa):
            axes=fig.axes
            axes[0].scatter(actual_T, actual_afa, c='black', zorder=2)
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def clearfig():
            axes=fig.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_agg):
            figure_canvas_agg.get_tk_widget().forget()
            plt.close('all')

        while True:
            event3, values3 = window3.read()

            if event3 == None:
                break

            if event3 == '-confirm1-':
                try:
                    liedu = values3['-combo1-']
                    shefang = values3['-combo2-']
                    fenzu = values3['-combo3-']
                    changdi = values3['-combo4-']
                    damping = float(values3['-damping-'])
                    
                    if liedu=='6度' and shefang=='多遇地震':
                        afa_max = 0.04
                    elif liedu=='6度' and shefang=='设防地震':
                        afa_max = 0.12
                    elif liedu=='6度' and shefang=='罕遇地震':
                        afa_max = 0.28
                    elif liedu=='6度' and shefang=='极罕遇地震':
                        afa_max = 0.36
                    elif liedu=='7度' and shefang=='多遇地震':
                        afa_max = 0.08
                    elif liedu=='7度' and shefang=='设防地震':
                        afa_max = 0.23
                    elif liedu=='7度' and shefang=='罕遇地震':
                        afa_max = 0.50
                    elif liedu=='7度' and shefang=='极罕遇地震':
                        afa_max = 0.7
                    elif liedu=='7度半' and shefang=='多遇地震':
                        afa_max = 0.12
                    elif liedu=='7度半' and shefang=='设防地震':
                        afa_max = 0.34
                    elif liedu=='7度半' and shefang=='罕遇地震':
                        afa_max = 0.72
                    elif liedu=='7度半' and shefang=='极罕遇地震':
                        afa_max = 1.0
                    elif liedu=='8度' and shefang=='多遇地震':
                        afa_max = 0.16
                    elif liedu=='8度' and shefang=='设防地震':
                        afa_max = 0.45
                    elif liedu=='8度' and shefang=='罕遇地震':
                        afa_max = 0.90
                    elif liedu=='8度' and shefang=='极罕遇地震':
                        afa_max = 1.35
                    elif liedu=='8度半' and shefang=='多遇地震':
                        afa_max = 0.24
                    elif liedu=='8度半' and shefang=='设防地震':
                        afa_max = 0.68
                    elif liedu=='8度半' and shefang=='罕遇地震':
                        afa_max = 1.2
                    elif liedu=='8度半' and shefang=='极罕遇地震':
                        afa_max = 2.0
                    elif liedu=='9度' and shefang=='多遇地震':
                        afa_max = 0.32
                    elif liedu=='9度' and shefang=='设防地震':
                        afa_max = 0.90
                    elif liedu=='9度' and shefang=='罕遇地震':
                        afa_max = 1.4
                    else:
                        afa_max = 2.7


                    if fenzu=='第一组' and changdi=='I0':
                        tg = 0.20
                    elif fenzu=='第一组' and changdi=='I1':
                        tg = 0.25
                    elif fenzu=='第一组' and changdi=='Ⅱ':
                        tg = 0.35
                    elif fenzu=='第一组' and changdi=='Ⅲ':
                        tg = 0.45
                    elif fenzu=='第一组' and changdi=='Ⅳ':
                        tg = 0.65
                    elif fenzu=='第二组' and changdi=='I0':
                        tg = 0.25
                    elif fenzu=='第二组' and changdi=='I1':
                        tg = 0.30
                    elif fenzu=='第二组' and changdi=='Ⅱ':
                        tg = 0.40
                    elif fenzu=='第二组' and changdi=='Ⅲ':
                        tg = 0.55
                    elif fenzu=='第二组' and changdi=='Ⅳ':
                        tg = 0.75
                    elif fenzu=='第三组' and changdi=='I0':
                        tg = 0.30
                    elif fenzu=='第三组' and changdi=='I1':
                        tg = 0.35
                    elif fenzu=='第三组' and changdi=='Ⅱ':
                        tg = 0.45
                    elif fenzu=='第三组' and changdi=='Ⅲ':
                        tg = 0.65
                    elif fenzu=='第三组' and changdi=='Ⅳ':
                        tg = 0.90

                    gama = 0.9 + (0.05-damping)/(0.3+6*damping)
                    eta1 = 0.02 + (0.05-damping)/(4+32*damping)
                    if eta1<0:
                        eta1=0
                    else:
                        eta1=eta1
                    eta2 = 1 + (0.05-damping)/(0.08+1.6*damping)
                    if eta2<0.55:
                        eta2=0.55
                    else:
                        eta2=eta2

                    times, afas = figdata(afa_max, tg, gama, eta1, eta2)
                    clearfig()
                    updatefig(times, afas, 1, label='设计反应谱')
                except:
                    sg.popup('请输入结构阻尼比', title='提示', text_color='red', font=20)
            
            if event3=='-delete-':
                window3['-combo1-'].update('6度')
                window3['-combo2-'].update('设防地震')
                window3['-combo3-'].update('第一组')
                window3['-combo4-'].update('Ⅱ')
                window3['-damping-'].update('0.05')
                window3['-period-'].update('')
                window3['-sa-'].update('')
                delete_fig_agg(figure_canvas_agg)
                sg.popup('清除成功！', font=20)

            if event3=='-confirm2-':
                try:
                    actual_T = float(values3['-period-'])
                    for id in range(len(times)-1):
                        if actual_T>times[id] and actual_T<=times[id+1]:
                            actual_afa = (afas[id+1]-afas[id])/(times[id+1]-times[id])*(actual_T-times[id])+afas[id]
                    window3['-sa-'].update(f'{round(actual_afa, 5)}')
                    clearfig()
                    updatefig(times, afas, 1, label='设计反应谱')
                    scatterfig(actual_T, actual_afa)
                except:
                    sg.popup("请输入正确的周期(0至6s)", title='提示', text_color='red', font=20)
            
            if event3=='-save-':
                try:
                    file_location = values3['-savelocation-']
                    with open(file_location, 'w') as f:
                        f.writelines(f'地震烈度: {liedu}\n')
                        f.writelines(f'设防等级: {shefang}\n')
                        f.writelines(f'地震分组: {fenzu}\n')
                        f.writelines(f'场地类别: {changdi}\n')
                        f.writelines(f'结构阻尼比: {damping}\n')
                        f.writelines(f'水平地震影响系数: {afa_max}\n')
                        f.writelines(f'\n')
                        f.writelines(f'周期 谱加速度\n')
                        for time, afa in zip(times, afas):
                            f.writelines("%f %f\n" %(time, afa))
                    sg.popup('保存成功！', title='提示', font=20)
                except:
                    sg.popup("请先选择要保存的文件名!", title='提示', text_color='red', font=20)
            
            if event3 == sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
                    window3.close()
                    window0.close()
    
            if event3 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提醒', text_color='red', font=20)=="Yes":
                    window3.close()
                    win3_active=False
                    window0.UnHide()
            

                
# ## -----------window4----------window4----------window4----------window4----------window4----------window4----------window4----------                
    if event0 == '-earthquakewaveprocessing-' and not win4_active:
        win4_active=True
        window0.Hide()
        formxz1 = ['PEER地震波', 'Column_2', 'Column_1']
        formxz2 = ['Column_2', 'Column_1']
        morenxuanze1 = ['按默认', '自定义']
        morenxuanze2 =['按默认', '调整PGA至目标值']

        table_content = []
        total_filenames_location = []

        layout41 = [
            [sg.Text('导入地震波格式', font=20), sg.Combo(formxz1, key='-combo1-', size=(12, 1), font=20, enable_events=True, default_value='PEER地震波'), sg.Text("导入地震波时间间隔", key='-text1-', visible=False, font=20), sg.In("", key="-deltatime1-", size=(12, 1), visible=False, font=20)],
            [sg.FilesBrowse("添加地震波文件", font=20, file_types=(('ALL Files', '*.txt'), ('ALL Files', '*.AT2')), target='-filenames-'), sg.In("", key='-filenames-', font=20), sg.B('添加', key='-confirm1-', font=20), sg.B("清除选择", key="-clearchoose-", font=20), sg.B("删除", key="-delete-", font=20), sg.B("绘图", key="-draw-", font=20)],
            [sg.Table(headings=['序号', '地震名称', '时长', 'PGA', "导入格式", '时间间隔'],
                    values=table_content,
                    expand_x=True,
                    expand_y=True,
                    key='-table-', 
                    font=20,
                    justification='center',
                    select_mode="extended",
                    enable_events=True)],
            [sg.Canvas(key='-toolbar-', size=(80, 1))],
            [sg.Canvas(key='-canvas-', expand_x=True, expand_y=True)]
        ]

        layout42 = [
            [sg.Text("批量导出地震波格式", font=20), sg.Combo(formxz2, key='-combo2-', font=20, default_value='Column_2')],
            [sg.Text('批量导出地震波时间间隔', font=20), sg.Combo(morenxuanze1, key='-combo3-', font=20, size=(8, 1), enable_events=True, default_value='按默认'), sg.Text("导出地震波时间间隔", key='-text2-', visible=False, font=20), sg.In('', key='-deltatime2-', visible=False, size=(8, 1), font=20)], 
            [sg.Text('批量调整峰值加速度', font=20), sg.Combo(morenxuanze2, key='-combo4-', font=20, size=(15, 1), enable_events=True, default_value='按默认'), sg.Text("目标PGA", key='-text3-', visible=False, font=20), sg.In('', key='-aimpga-', visible=False, size=(8, 1), font=20)],
            [sg.FolderBrowse('选择保存文件夹', target='-savelocation-', font=20), sg.In('', key='-savelocation-', font=20), sg.B('保存', key='-save-', font=20), sg.B('清除所有', key='-clearall-', font=20), sg.B('返回至初始界面', key='-back-', button_color=('red', 'skyblue'), font=20)]
        ]

        layout40 = [
            [sg.Frame('地震波选择', layout41, title_color='red', font=('楷体', 15), expand_x=True, expand_y=True, element_justification='center')],
            [sg.Frame('峰值加速度调幅与批量导出', layout42, title_color='red', font=('楷体', 15), expand_x=True, element_justification='center', expand_y=True)],
        ]

        layout4 = [[sg.Column(layout40, scrollable=True, vertical_scroll_only=True, size=(1200, 1100))]]

        window4 = sg.Window("土木工程数据常用处理软件 ----- 地震波格式转换&峰值调整", layout4, location=(500, 20), finalize=True, resizable=True, enable_close_attempted_event=True)


        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)
        fig = matplotlib.figure.Figure()
        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False
        fig.add_subplot(111)
        axes=fig.axes
        axes[0].set_xlabel('时间')
        axes[0].set_ylabel('加速度')
        figure_canvas_agg = FigureCanvasTkAgg(fig, window4['-canvas-'].TKCanvas)
        toolbar = Toolbar(figure_canvas_agg, window4['-toolbar-'].TKCanvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack()

        def updatefig(times, accs, alpha, label):
            axes=fig.axes
            axes[0].plot(times, accs, alpha=alpha, label=label)
            axes[0].set_xlim(0, )
            axes[0].set_xlabel('时间')
            axes[0].set_ylabel('加速度')
            axes[0].legend()
            toolbar.update()
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def clearfig():
            axes=fig.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_agg):
            figure_canvas_agg.get_tk_widget().forget()
            plt.close('all')


        while True:
            event4, values4 = window4.read()

            if event4 == None:
                break
            
            if event4 == '-combo1-': 
                if values4['-combo1-'] == 'Column_1':
                    window4['-text1-'].update(visible=True)
                    window4['-deltatime1-'].update('', visible=True)
                else:
                    window4['-text1-'].update(visible=False)
                    window4['-deltatime1-'].update('', visible=False)

            if event4 == '-confirm1-':
                try:
                    window4['-filenames-'].update('')
                    filenames = str(values4['-filenames-']).split(';')
                    current_filenames = []
                    for i in range(len(table_content)):
                        current_filenames.append(table_content[i][1])
                    for filename in filenames:
                        earthquakename = filename.split('/')[-1]

                        if earthquakename in current_filenames:
                            sg.popup('请勿重复添加地震波', title='提示', text_color='red', font=20)
                            break
                        
                        else:
                            total_filenames_location.append(filename)

                            if values4['-combo1-'] == 'PEER地震波':
                                with open(filename, 'r') as f:
                                    lines = f.readlines()
                                strs = lines[3]
                                deltatime_input = float(strs.split('=')[2].split(' ')[3])
                                datas = pd.read_table(filename, header=None, skiprows=4, sep='\t')
                                accs = []
                                for i in range(len(datas)):
                                    col = datas.values[i]
                                    zifu = str(col).split('\'')[1]
                                    for each_v in zifu.split(' '):
                                        try:
                                            accs.append(float(each_v))
                                        except:
                                            continue
                
                                total_num = len(accs)
                                time = (total_num-1)*deltatime_input
                                accs = np.array(accs)
                                PGAvalue = max(abs(accs))
                                table_content.append([len(table_content)+1, earthquakename, time, PGAvalue, values4["-combo1-"], deltatime_input])
                                window4['-table-'].update(table_content)

                            if values4['-combo1-'] == 'Column_2':
                                datas = np.loadtxt(filename, dtype=float, unpack=True)
                                deltatime_input = datas[0][1]-datas[0][0]
                                time = datas[0][-1]
                                PGAvalue = max(abs(datas[1]))
                                table_content.append([len(table_content)+1, earthquakename, time, PGAvalue, values4["-combo1-"], deltatime_input])
                                window4['-table-'].update(table_content)

                            if values4['-combo1-'] == 'Column_1':
                                try:
                                    deltatime_input = float(values4['-deltatime1-'])
                                    datas = np.loadtxt(filename, dtype=float, unpack=True)
                                    time = (len(datas)-1)*deltatime_input
                                    PGAvalue = max(abs(datas))
                                    table_content.append([len(table_content)+1, earthquakename, time, PGAvalue, values4["-combo1-"], deltatime_input])
                                    window4['-table-'].update(table_content)
                                except:
                                    sg.popup('请先输入地震波时间间隔', title='提示', text_color='red', font=20)
                
                except:
                    sg.popup('请先检查输入地震波格式是否一致', title='提示', text_color='red', font=20)
                
            if event4 == "-clearchoose-":
                window4["-filenames-"].update("")
            
            if event4 == "-delete-":
                index = values4["-table-"]
                num = -1
                for i in index:
                    num += 1
                    id = i - num
                    table_content.pop(id)
                    total_filenames_location.pop(id)
                for id in range(len(table_content)):
                    table_content[id][0] = id+1
                window4['-table-'].update(table_content)

            if event4 == "-draw-":
                clearfig()
                index=values4["-table-"]
                
                if len(index)==0:
                    sg.popup("请先在下表中点击需要绘图的地震波", title="提示", text_color="red", font=20)
                
                if len(index)>1:
                    sg.popup("每次仅限绘制一条地震波时程曲线", title="提示", text_color="red", font=20)

                if len(index)==1:
                    id=index[0]
                    type_wave = table_content[id][4]
                    choose_filename_location = total_filenames_location[id]
                    
                    if type_wave=='PEER地震波':
                        deltatime_input = float(table_content[id][5])
                        datas = pd.read_table(choose_filename_location, header=None, skiprows=4, sep='\t')
                        accs = []
                        for i in range(len(datas)):
                            col = datas.values[i]
                            zifu = str(col).split('\'')[1]
                            for each_v in zifu.split(' '):
                                try:
                                    accs.append(float(each_v))
                                except:
                                    continue
                        total_num = len(accs)
                        time = (total_num-1)*deltatime_input
                        times = np.linspace(0, time, total_num)
                        accs = np.array(accs)
                        updatefig(times, accs, 1, label='原地震波')
                    
                    if type_wave=='Column_2':
                        datas = np.loadtxt(choose_filename_location, dtype=float, unpack=True)
                        times = datas[0]
                        accs = datas[1]
                        updatefig(times, accs, 1, label='原地震波')

                    if type_wave=='Column_1':
                        deltatime_input = float(table_content[id][5])
                        datas = np.loadtxt(choose_filename_location, dtype=float, unpack=True)
                        accs = datas
                        time = (len(datas)-1)*deltatime_input
                        times = np.linspace(0, time, len(datas))
                        updatefig(times, accs, 1, label='原地震波')

            if event4 == '-combo3-':
                if values4['-combo3-'] == '自定义':
                    window4['-text2-'].update(visible=True)
                    window4['-deltatime2-'].update('', visible=True)
                else:
                    window4['-text2-'].update(visible=False)
                    window4['-deltatime2-'].update('', visible=False)

            if event4 == '-combo4-':
                if values4['-combo4-'] == '调整PGA至目标值':
                    window4['-text3-'].update(visible=True)
                    window4['-aimpga-'].update('', visible=True)
                else:
                    window4['-text3-'].update(visible=False)
                    window4['-aimpga-'].update('', visible=False)

            if event4 == '-save-':
                    try:
                        savefolder_location = values4['-savelocation-']
                        for i in range(len(table_content)):
                            filename_location = total_filenames_location[i]
                            type_filewave = table_content[i][4]
                            deltatime_input = table_content[i][5]
                            total_time = table_content[i][2]
                            filename = table_content[i][1]
                            

                            if values4['-combo3-']=='自定义':
                                deltatime_output = float(values4['-deltatime2-'])
                            else:
                                deltatime_output = deltatime_input
                            
                            new_times, new_accs, new_accs1 = [], [], []
                            
                            if type_filewave=='PEER地震波':
                                datas = pd.read_table(filename_location, header=None, skiprows=4, sep='\t')
                                accs = []
                                for i in range(len(datas)):
                                    col = datas.values[i]
                                    zifu = str(col).split('\'')[1]
                                    for each_v in zifu.split(' '):
                                        try:
                                            accs.append(float(each_v))
                                        except:
                                            continue
                                total_num = len(accs)
                                times = np.linspace(0, total_time, total_num)
                                accs = np.array(accs)

                            if type_filewave=='Column_2':
                                datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                                times = datas[0]
                                accs = datas[1]

                            if type_filewave=='Column_1':
                                datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                                accs = datas
                                times = np.linspace(0, total_time, len(datas))
                            
                            total_new_num = int(total_time/deltatime_output)+1
                            new_times = np.linspace(0, deltatime_output*(total_new_num-1), total_new_num)

                            for ix in range(total_new_num):
                                for inx in range(len(times)-1):
                                    if times[inx]<=new_times[ix] and times[inx+1]>new_times[ix]:
                                        new_acc = (accs[inx+1]-accs[inx])/(times[inx+1]-times[inx])*(new_times[ix]-times[inx])+accs[inx]
                                        new_accs.append(new_acc)
                            
                            if total_time % deltatime_output<0.001:
                                new_accs.append(accs[-1])
                            else:
                                pass
                                        
                            PGA_original = float(np.max(np.abs(new_accs)))

                            if values4['-combo4-']=='调整PGA至目标值':
                                PGA_finally = float(values4['-aimpga-'])
                            else:
                                PGA_finally = PGA_original

                            
                            amp_factor = PGA_finally / PGA_original
                            for i in range(len(new_accs)):
                                new_accs1.append(amp_factor*new_accs[i])
                            savefile_location = savefolder_location+'/'+f'PGA为{PGA_finally}_'+filename
                            
                            with open(savefile_location, 'w') as f:
                                if values4['-combo2-'] == 'Column_2':
                                    for t, a in zip(new_times, new_accs1):
                                        f.writelines("%f %f\n" %(t, a))
                                
                                if values4['-combo2-'] == 'Column_1':
                                    for a in new_accs1:
                                        f.writelines('%f\n' %a)

                        sg.popup_timed("保存成功！", title="提示", font=20)
                    
                    except:
                        sg.popup("请先选择目标保存文件夹!", title="提示", text_color='red', font=20)
            
            if event4 == '-clearall-':
                window4['-combo1-'].update('PEER地震波')
                window4['-text1-'].update(visible=False)
                window4['-deltatime1-'].update('', visible=False)
                window4['-filenames-'].update('')
                table_content=[]
                window4['-table-'].update(table_content)
                delete_fig_agg(figure_canvas_agg)
                total_filenames_location=[]
                window4['-combo2-'].update('Column_2')
                window4['-combo3-'].update('按默认')
                window4['-deltatime2-'].update('', visible=False)
                window4['-text2-'].update(visible=False)
                window4['-text3-'].update(visible=False)
                window4['-aimpga-'].update('', visible=False)
                window4['-savelocation-'].update('')
                sg.popup('清除成功！', font=20)

            if event4 == sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
                    window4.close()
                    window0.close()
    
            if event4 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提醒', text_color='red', font=20)=="Yes":
                    window4.close()
                    win4_active=False
                    window0.UnHide()
     


# ## -----------window5----------window5----------window5----------window5----------window5----------window5----------window5----------    
    if event0 == '-earthquakewavefilter-' and not win5_active:
        win5_active=True
        window0.Hide()

        table_content = []
        total_filenames_location = []
        needvd = ["需要", '不需要']

        layout511 = [
            [sg.Text('原地震波加速度时域曲线', font=20), sg.Text('                                                           '), sg.Text('原地震波加速度频域曲线', font=20), sg.Text('                                                           '), sg.Text('原地震波位移时域曲线', font=20)],
            [sg.Canvas(key='-toolbar1-', expand_x=True, expand_y=True), sg.Canvas(key='-toolbar2-', expand_x=True, expand_y=True), sg.Canvas(key='-toolbar5-', expand_x=True, expand_y=True)],
            [sg.Canvas(key='-canvas1-', expand_x=True, expand_y=True), sg.Canvas(key='-canvas2-', expand_x=True, expand_y=True), sg.Canvas(key='-canvas5-', expand_x=True, expand_y=True)],
            [sg.Text('滤波后地震波加速度时域曲线', font=20), sg.Text('                                                '), sg.Text('滤波后地震波加速度频域曲线', font=20), sg.Text('                                                '), sg.Text('滤波后地震波位移时域曲线', font=20)],
            [sg.Canvas(key='-toolbar3-', expand_x=True, expand_y=True), sg.Canvas(key='-toolbar4-', expand_x=True, expand_y=True), sg.Canvas(key='-toolbar6-', expand_x=True, expand_y=True)],
            [sg.Canvas(key='-canvas3-', expand_x=True, expand_y=True), sg.Canvas(key='-canvas4-', expand_x=True, expand_y=True), sg.Canvas(key='-canvas6-', expand_x=True, expand_y=True)],
        ]

        layout51 = [
            [sg.FilesBrowse("添加地震波文件", key='-openfiles-', font=20, file_types=(('ALL Files', '*.txt'), ('ALL Files', '*.AT2')), target='-filenames-'), sg.In("", key='-filenames-', font=20), sg.B('添加', key='-confirm1-', font=20), sg.B("清除选择", key="-clearchoose-", font=20), sg.B("删除", key="-delete-", font=20)],
            [sg.Text('截止频率', font=20), sg.In('nan', size=(12, 1), font=20, key='-cutfrequency-'), sg.B('确定并绘图', key='-confirm2-', font=20)],
            [sg.Table(headings=['序号', '地震名称', '时长', 'PGA', '截止频率'],
                    values=table_content,
                    expand_x=True,
                    expand_y=True, 
                    font=20,
                    num_rows=5,
                    key='-table-',
                    justification='center',
                    select_mode="extended",
                    enable_events=True)],
            [sg.Frame('地震波时程记录', layout511, title_color='black', font=('楷体'), expand_x=True, expand_y=True, element_justification='center')]
        ]

        layout52 = [
            [sg.Text("是否需要导出滤波后的地面运动速度和位移数据", font=20), sg.Combo(needvd, default_value='不需要', size=(12, 1), font=20, key='-combo1-')],
            [sg.FolderBrowse('选择保存文件夹', target='-savelocation-', font=20), sg.In('', key='-savelocation-', visible=True, font=20), sg.B('保存', key='-save-', font=20), sg.B('清除所有', key='-clearall-', font=20), sg.B('返回至初始界面', key='-back-', button_color=('red', 'skyblue'), font=20)]
        ]

        layout50 = [
            [sg.Frame('地震波选择与滤波处理', layout51, title_color='red', font=('楷体'), expand_x=True, expand_y=True, element_justification='center')],
            [sg.Frame('批量导出加速度、速度与位移时程数据', layout52, title_color='red', font=('楷体'), expand_x=True, element_justification='center', expand_y=True)]
        ]

        layout5 = [[sg.Column(layout50, scrollable=True, vertical_scroll_only=True, size=(1500, 1200))]]

        window5 = sg.Window("土木工程数据常用处理软件 ----- 地震动滤波处理&时程转换", layout5, size=(1500, 1200), location=(350, 0), finalize=True, resizable=True, enable_close_attempted_event=True)

        def accs_fft(times, accs):
            freqs = nf.fftfreq(times.size, times[1]-times[0])
            complex_accs = nf.fft(accs)
            pows = np.abs(complex_accs)
            return freqs, pows, complex_accs
        
        def accs_fft_cut(freqs, complex_accs, cut_frequency):
            complex_filter = complex_accs.copy()
            complex_filter[abs(freqs)<cut_frequency] = 0
            pows_filter = np.abs(complex_filter)
            return freqs, pows_filter, complex_filter
        
        def accs_ifft(complex_filter):
            accs_filter = np.fft.ifft(complex_filter).real
            return accs_filter
        
        def accum_velocity(times, accs):
            vels = []
            for id in range(len(times)):
                v = integrate.trapz(accs[:id], times[:id])
                vels.append(v)
            vels =np.array(vels)
            return vels
        
        def accum_displacement(times, vels):
            disps = []
            for id in range(len(times)):
                d = integrate.trapz(vels[:id], times[:id])
                disps.append(d)
            disps =np.array(disps)
            return disps

        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)

        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False

        fig1 = matplotlib.figure.Figure(figsize=(4.7, 2.8))
        fig1.add_subplot(111)
        axes1 = fig1.axes
        axes1[0].set_xlabel('时间')
        axes1[0].set_ylabel('加速度')

        fig2 = matplotlib.figure.Figure(figsize=(4.7, 2.8))
        fig2.add_subplot(111)
        axes2 = fig2.axes
        axes2[0].set_xlabel('幅值')
        axes2[0].set_ylabel('频率')

        fig3 = matplotlib.figure.Figure(figsize=(4.7, 2.8))
        fig3.add_subplot(111)
        axes3 = fig3.axes
        axes3[0].set_xlabel('时间')
        axes3[0].set_ylabel('加速度')

        fig4 = matplotlib.figure.Figure(figsize=(4.7, 2.8))
        fig4.add_subplot(111)
        axes4 = fig2.axes
        axes4[0].set_xlabel('幅值')
        axes4[0].set_ylabel('频率')

        fig5 = matplotlib.figure.Figure(figsize=(4.7, 2.8))
        fig5.add_subplot(111)
        axes5 = fig5.axes
        axes5[0].set_xlabel('时间')
        axes5[0].set_ylabel('位移')

        fig6 = matplotlib.figure.Figure(figsize=(4.7, 2.8))
        fig6.add_subplot(111)
        axes6 = fig6.axes
        axes6[0].set_xlabel('时间')
        axes6[0].set_ylabel('位移')
        figure_canvas_agg1 = FigureCanvasTkAgg(fig1, window5['-canvas1-'].TKCanvas)
        figure_canvas_agg2 = FigureCanvasTkAgg(fig2, window5['-canvas2-'].TKCanvas)
        figure_canvas_agg3 = FigureCanvasTkAgg(fig3, window5['-canvas3-'].TKCanvas)
        figure_canvas_agg4 = FigureCanvasTkAgg(fig4, window5['-canvas4-'].TKCanvas)
        figure_canvas_agg5 = FigureCanvasTkAgg(fig5, window5['-canvas5-'].TKCanvas)
        figure_canvas_agg6 = FigureCanvasTkAgg(fig6, window5['-canvas6-'].TKCanvas)
        toolbar1 = Toolbar(figure_canvas_agg1, window5['-toolbar1-'].TKCanvas)
        toolbar2 = Toolbar(figure_canvas_agg2, window5['-toolbar2-'].TKCanvas)
        toolbar5 = Toolbar(figure_canvas_agg5, window5['-toolbar5-'].TKCanvas)
        toolbar3 = Toolbar(figure_canvas_agg3, window5['-toolbar3-'].TKCanvas)
        toolbar4 = Toolbar(figure_canvas_agg4, window5['-toolbar4-'].TKCanvas)
        toolbar6 = Toolbar(figure_canvas_agg6, window5['-toolbar6-'].TKCanvas)
        
        figure_canvas_agg1.draw()
        figure_canvas_agg1.get_tk_widget().pack()

        figure_canvas_agg2.draw()
        figure_canvas_agg2.get_tk_widget().pack()
        
        figure_canvas_agg3.draw()
        figure_canvas_agg3.get_tk_widget().pack()
        
        figure_canvas_agg4.draw()
        figure_canvas_agg4.get_tk_widget().pack()
        
        figure_canvas_agg5.draw()
        figure_canvas_agg5.get_tk_widget().pack()
        
        figure_canvas_agg6.draw()
        figure_canvas_agg6.get_tk_widget().pack()

        def updatefig_time(figname, figure_canvas_aggname, toolbar, times, accs, alpha, label):
            axes=figname.axes
            axes[0].plot(times, accs, alpha=alpha, label=label)
            axes[0].legend()
            axes[0].set_xlim(0, )
            figure_canvas_aggname.draw()
            toolbar.update()
            figure_canvas_aggname.get_tk_widget().pack()

        def updatefig_frequency(figname, figure_canvas_aggname, toolbar, frequencys, pows, alpha, label):
            axes=figname.axes
            axes[0].plot(frequencys, pows, alpha=alpha, label=label)
            axes[0].legend()
            axes[0].set_xlim(0, )
            figure_canvas_aggname.draw()
            toolbar.update()
            figure_canvas_aggname.get_tk_widget().pack()

        def clearfig(figname):
            axes=figname.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_aggname):
            figure_canvas_aggname.get_tk_widget().forget()
            plt.close('all')


        while True:
            event5, values5 = window5.read()

            if event5 == None:
                break

            if event5 == '-confirm1-':
                try:
                    window5['-filenames-'].update('')
                    filenames = str(values5['-filenames-']).split(';')
                    current_filenames = []
                    for i in range(len(table_content)):
                        current_filenames.append(table_content[i][1])
                    for filename in filenames:
                        earthquakename = filename.split('/')[-1]

                        if earthquakename in current_filenames:
                            sg.popup('请勿重复添加地震波', title='提示', text_color='red', font=20)
                            break

                        else:
                            total_filenames_location.append(filename)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            time = datas[0][-1]
                            PGAvalue = max(abs(datas[1]))
                            table_content.append([len(table_content)+1, earthquakename, time, PGAvalue, np.nan])
                            window5['-table-'].update(table_content)
                
                except:
                    sg.popup('请先选择要打开的地震波文件', title='提示', text_color='red', font=20)

            if event5 == '-clearchoose-':
                window5["-filenames-"].update("")

            if event5 == "-delete-":
                index = values5["-table-"]
                num = -1
                for i in index:
                    num += 1
                    id = i - num
                    table_content.pop(id)
                    total_filenames_location.pop(id)
                for id in range(len(table_content)):
                    table_content[id][0] = id+1
                window5['-table-'].update(table_content)

            if event5 == '-confirm2-':
                clearfig(fig1)
                clearfig(fig2)
                clearfig(fig3)
                clearfig(fig4)
                clearfig(fig5)
                clearfig(fig6)
                cut_frequency = values5['-cutfrequency-']
                index = values5['-table-']
                
                if len(index)==0:
                    sg.popup('请先在下表中点击一条需要处理的地震波', text_color='red', title='提示', font=20)
                
                elif len(index)>1:
                    sg.popup('每次仅支持处理一条地震波', text_color='red', title='提示', font=20)
                
                else:
                    filename = total_filenames_location[index[0]]
                    if cut_frequency == 'nan':
                        datas = np.loadtxt(filename, dtype=float, unpack=True)
                        times = datas[0]
                        accs = datas[1]
                        freqs, pows, complex_accs = accs_fft(times, accs)
                        vels = accum_velocity(times, accs)
                        disps = accum_displacement(times, vels)
                        updatefig_time(fig1, figure_canvas_agg1, toolbar1, times, accs, 1, label='原地震波')
                        updatefig_frequency(fig2, figure_canvas_agg2, toolbar2, freqs, pows, 1, label='原地震波')
                        updatefig_time(fig5, figure_canvas_agg5, toolbar5, times, disps, 1, label='原地震波')
                    
                    if cut_frequency != 'nan':
                        try:
                            cutfrequency = float(cut_frequency)
                            datas = np.loadtxt(filename, dtype=float, unpack=True)
                            times = datas[0]
                            accs = datas[1]
                            freqs, pows, complex_accs = accs_fft(times, accs)
                            vels = accum_velocity(times, accs)
                            disps = accum_displacement(times, vels)
                            freqs, pows_filter, complex_filter = accs_fft_cut(freqs, complex_accs, cutfrequency)
                            accs_filter = accs_ifft(complex_filter)
                            vels_filter = accum_velocity(times, accs_filter)
                            disps_filter = accum_displacement(times, vels_filter)
                            updatefig_time(fig1, figure_canvas_agg1, toolbar1, times, accs, 1, label='原地震波')
                            updatefig_frequency(fig2, figure_canvas_agg2, toolbar2, freqs, pows, 1, label='原地震波')
                            updatefig_frequency(fig4, figure_canvas_agg4, toolbar4, freqs, pows_filter, 1, label='滤波后')
                            updatefig_time(fig3, figure_canvas_agg3, toolbar3, times, accs_filter, 1, label='滤波后')
                            updatefig_time(fig5, figure_canvas_agg5, toolbar5, times, disps, 1, label='原地震波')
                            updatefig_time(fig6, figure_canvas_agg6, toolbar6, times, disps_filter, 1, label='滤波后')
                            table_content[index[0]][4] = cutfrequency
                            window5['-table-'].update(table_content)
                        
                        except:
                            sg.popup("请输入正确的截止频率", title='提示', text_color='red', font=20)
                            
            if event5 == '-save-':
                savefolder_location = values5['-savelocation-']

                if values5['-combo1-'] == '不需要':
                    try:
                        for i in range(len(table_content)):
                                filename_location = total_filenames_location[i]
                                filename = table_content[i][1]
                                cut_frequency_ = table_content[i][4]

                                if math.isnan(cut_frequency_):
                                    savefile_location = savefolder_location+f'/{i}.未滤波加速度_'+filename
                                    datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                                    times = datas[0]
                                    accs = datas[1]
                                    with open(savefile_location, 'w') as f:
                                        for t, a in zip(times, accs):
                                            f.writelines('%f %f\n' %(t, a))
                                else:
                                    savefile_location = savefolder_location+f'/{i}.滤波后加速度_'+filename
                                    datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                                    times = datas[0]
                                    accs = datas[1]
                                    freqs, pows, complex_accs = accs_fft(times, accs)
                                    freqs, pows_filter, complex_filter = accs_fft_cut(freqs, complex_accs, cutfrequency)
                                    accs_filter = accs_ifft(complex_filter)
                                    with open(savefile_location, 'w') as f:
                                        for t, a in zip(times, accs_filter):
                                            f.writelines('%f %f\n' %(t, a))
                        sg.popup('保存成功', title='提示', text_color='green', font=20)
                    except:
                        sg.popup('请先选择保存的文件夹！', title='提示', text_color='red', font=20)

                if values5['-combo1-'] == '需要':
                    try:
                        for i in range(len(table_content)):
                                filename_location = total_filenames_location[i]
                                filename = table_content[i][1]
                                cut_frequency_ = table_content[i][4]

                                if math.isnan(cut_frequency_):
                                    savefile_location_a = savefolder_location+f'/{i}.未滤波加速度_'+filename
                                    savefile_location_v = savefolder_location+f'/{i}.未滤波速度_'+filename
                                    savefile_location_d = savefolder_location+f'/{i}.未滤波位移_'+filename
                                    datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                                    times = datas[0]
                                    accs = datas[1]
                                    vels = accum_velocity(times, accs)
                                    disps = accum_displacement(times, vels)

                                    with open(savefile_location_a, 'w') as f:
                                        for t, a in zip(times, accs):
                                            f.writelines('%f %f\n' %(t, a))
                                    
                                    with open(savefile_location_v, 'w') as f:
                                        for t, v in zip(times, vels):
                                            f.writelines('%f %f\n' %(t, v))

                                    with open(savefile_location_d, 'w') as f:
                                        for t, d in zip(times, disps):
                                            f.writelines('%f %f\n' %(t, d))
                                
                                else:
                                    savefile_location_a = savefolder_location+f'/{i}.滤波后加速度_'+filename
                                    savefile_location_v = savefolder_location+f'/{i}.滤波后速度_'+filename
                                    savefile_location_d = savefolder_location+f'/{i}.滤波后位移_'+filename
                                    datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                                    times = datas[0]
                                    accs = datas[1]
                                    vels = accum_velocity(times, accs)
                                    disps = accum_displacement(times, vels)
                                    freqs, pows, complex_accs = accs_fft(times, accs)
                                    freqs, pows_filter, complex_filter = accs_fft_cut(freqs, complex_accs, cutfrequency)
                                    accs_filter = accs_ifft(complex_filter)
                                    vels_filter = accum_velocity(times, accs_filter)
                                    disps_filter = accum_displacement(times, vels_filter)
                                    
                                    with open(savefile_location_a, 'w') as f:
                                        for t, a in zip(times, accs_filter):
                                            f.writelines('%f %f\n' %(t, a))
                                    
                                    with open(savefile_location_v, 'w') as f:
                                        for t, v in zip(times, vels_filter):
                                            f.writelines('%f %f\n' %(t, v))

                                    with open(savefile_location_d, 'w') as f:
                                        for t, d in zip(times, disps_filter):
                                            f.writelines('%f %f\n' %(t, d))
                        sg.popup('保存成功', title='提示', text_color='green', font=20)
                    except:
                        sg.popup('请先选择保存的文件夹！', title='提示', text_color='red', font=20)
                    

            if event5 == '-clearall-':
                window5['-filenames-'].update('')
                window5['-cutfrequency-'].update('nan')
                table_content = []
                window5['-table-'].update(table_content)
                delete_fig_agg(figure_canvas_agg1)
                delete_fig_agg(figure_canvas_agg2)
                delete_fig_agg(figure_canvas_agg3)
                delete_fig_agg(figure_canvas_agg4)
                delete_fig_agg(figure_canvas_agg5)
                delete_fig_agg(figure_canvas_agg6)
                total_filenames_location = []
                window5['-combo1-'].update('不需要')
                window5['-savelocation-'].update('')
                sg.popup('清除成功！', title='提示', font=20)

            if event5 == sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提示', text_color='red', font=20)=="Yes":
                    window5.close()
                    window0.close()

            if event5 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提示', text_color='red', font=20)=="Yes":
                    window5.close()
                    win5_active=False
                    window0.UnHide()


                            
# ## -----------window6----------window6----------window6----------window6----------window6----------window6----------window6----------  
    if event0 == '-fragility-' and not win6_active:
        win6_active = True
        window0.Hide()

        table_content = []
        limitstatenum = ["1", "2", "3", "4", "5", "6"]
        nihenum = ["一次", '二次']
        
        layout62 = [
            [sg.Text('IDA曲线与回归拟合公式', font=20), sg.T('                                                                                                                      '), sg.Text('易损性曲线绘制', font=20)],
            [sg.Combo(nihenum, key='-nihenum-', font=20, enable_events=True, size=(4, 1), default_value='一次'), sg.Text('多项式拟合', font=20), sg.Input('ln(DM) = b + a·ln(IM)', key='-equal1-', font=20, size=(44, 1), disabled=True), sg.Text('易损性曲线横坐标IM取值范围', font=20), sg.Input('', key='-startcoord-', font=20, size=(5, 1)), sg.Text('至', font=20), sg.Input('', key='-endcoord-', font=20, size=(5, 1))],
            [sg.Canvas(key='-toolbar1-', expand_x=True, expand_y=True), sg.Canvas(key='-toolbar2-', expand_x=True, expand_y=True)],
            [sg.Canvas(key='-canvas1-', expand_x=True, expand_y=True), sg.Canvas(key='-canvas2-', expand_x=True, expand_y=True)],
            [sg.Text('超越概率查询:     地震动强度IM', font=20), sg.In('', key='-imquery-', font=20, size=(8, 1)), sg.B('查询', key='-query-', font=20)],
            [sg.Text('超越概率:  LS1', font=20), sg.In('', key='-ls1query-', font=20, size=(8, 1), disabled=True), sg.Text('LS2', font=20, key='-textquery2-', visible=True), sg.In('', key='-ls2query-', font=20, size=(8, 1), visible=True, disabled=True), sg.Text('LS3', font=20, key='-textquery3-', visible=True), sg.In('', key='-ls3query-', font=20, size=(8, 1), visible=True, disabled=True), sg.Text('LS4', font=20, key='-textquery4-', visible=True), sg.In('', key='-ls4query-', font=20, size=(8, 1), visible=True, disabled=True), sg.Text('LS5', font=20, key='-textquery5-', visible=False), sg.In('', key='-ls5query-', font=20, size=(8, 1), visible=False, disabled=True), sg.Text('LS6', font=20, key='-textquery6-', visible=False), sg.In('', key='-ls6query-', font=20, size=(8, 1), visible=False, disabled=True)]
        ]

        layout61 = [
            [sg.Text('地震动强度指标IM', font=20), sg.In('', key='-im-', font=20, size=(10, 1)), sg.Text('结构需求指标DM', font=20), sg.In('', key='-dm-', font=20, size=(10, 1)), sg.B('添加', key='-confirm1-', font=20), sg.B('删除', key='-delete1-', font=20)],
            [sg.Table(headings=['序号', 'IM', 'DM'],
                    values=table_content,
                    expand_x=True,
                    expand_y=True, 
                    font=20,
                    num_rows=5,
                    key='-table-',
                    justification='center',
                    select_mode="extended",
                    enable_events=True)],
            [sg.Text('请选择损伤状态的数量', font=20), sg.Combo(limitstatenum, key='-limitstatenum-', enable_events=True, default_value='4', font=20, size=(10, 1))],
            [sg.Text('各损伤状态对应结构地震需求限值', font=20), sg.B('确定并绘图', key='-confirm2-', font=20)],
            [sg.Text('LS1限值', key='-textls1-', font=20), sg.In('', key='-ls1-', font=20, visible=True), sg.Text('LS2限值', key='-textls2-', font=20), sg.In('', key='-ls2-', font=20, visible=True)],
            [sg.Text('LS3限值', key='-textls3-', font=20), sg.In('', key='-ls3-', font=20, visible=True), sg.Text('LS4限值', key='-textls4-', font=20), sg.In('', key='-ls4-', font=20, visible=True)],
            [sg.Text('LS5限值', key='-textls5-', visible=False, font=20), sg.In('', key='-ls5-', font=20, visible=False), sg.Text('LS6限值', key='-textls6-', font=20, visible=False), sg.In('', key='-ls6-', font=20, visible=False)]
        ]
        

        layout60 = [
            [sg.Frame('数据填入', layout61, title_color='red', font=('楷体', 15), expand_x=True, expand_y=True)],
            [sg.Frame('易损性曲线', layout62, title_color='red', font=('楷体', 15), expand_x=True, expand_y=True)],
            [sg.FileSaveAs('选择保存文件路径', target='-savelocation-', font=20, file_types=(('*.txt', '*.txt'),)), sg.In('', key='-savelocation-', font=20), sg.B('保存', key='-save-', font=20), sg.B('清除所有', key='-clearall-', font=20), sg.B('返回至初始界面', key='-back-', button_color=('red', 'skyblue'), font=20)]
        ]
        layout6 = [[sg.Column(layout60, scrollable=True, vertical_scroll_only=True, size=(1400, 1100))]]

        window6 = sg.Window('土木工程数据常用处理软件 ----- IDA方法易损性分析', layout6, location=(400, 20), finalize=True, resizable=True, enable_close_attempted_event=True)


        def fragility_ida(xs, ys, deg):
            if deg == 2:
                poly = np.polyfit(xs, ys, deg=deg)
                a = poly[0]
                b = poly[1]
                c = poly[2]
                return a, b, c
            if deg == 1:
                poly = np.polyfit(xs, ys, deg=deg)
                a = poly[0]
                b = poly[1]
                return a, b

        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)

        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False
        fig1 = matplotlib.figure.Figure(figsize=(6.5, 4.3))
        fig1.add_subplot(111)
        axes1=fig1.axes
        axes1[0].set_xlabel('ln(IM)')
        axes1[0].set_ylabel('ln(DM)')
        figure_canvas_agg1 = FigureCanvasTkAgg(fig1, window6['-canvas1-'].TKCanvas)
        toolbar1 = Toolbar(figure_canvas_agg1, window6['-toolbar1-'].TKCanvas)
        figure_canvas_agg1.draw()
        figure_canvas_agg1.get_tk_widget().pack()

        fig2 = matplotlib.figure.Figure(figsize=(6.5, 4.3))
        fig2.add_subplot(111)
        axes2=fig2.axes
        axes2[0].set_xlabel('IM')
        axes2[0].set_ylabel('超越概率Pf/%')
        figure_canvas_agg2 = FigureCanvasTkAgg(fig2, window6['-canvas2-'].TKCanvas)
        toolbar2 = Toolbar(figure_canvas_agg2, window6['-toolbar2-'].TKCanvas)
        figure_canvas_agg2.draw()
        figure_canvas_agg2.get_tk_widget().pack()

        def updatefig(figname, figure_canvas_agg, toolbar, xs, ys, alpha, label):
            axes=figname.axes
            axes[0].plot(xs, ys, alpha=alpha, label=label, c='black')
            axes[0].legend()
            axes[0].set_xlabel("ln(IM)")
            axes[0].set_ylabel("ln(DM)")
            figure_canvas_agg.draw()
            toolbar.update()
            figure_canvas_agg.get_tk_widget().pack()

        def drawfig(figname, figure_canvas_agg, toolbar, xs, ys, alpha, label):
            axes=figname.axes
            col = (np.random.random(), np.random.random(), np.random.random())
            axes[0].plot(xs, ys, alpha=alpha, label=label, c=col)
            axes[0].legend()
            axes[0].set_xlabel('IM')
            axes[0].set_ylabel('超越概率Pf/%')
            figure_canvas_agg.draw()
            toolbar.update()
            figure_canvas_agg.get_tk_widget().pack()
        
        def scatterfig(figname, figure_canvas_agg, toolbar, xs, ys):
            axes=figname.axes
            axes[0].scatter(xs, ys, zorder=2, c=np.random.rand(len(xs),3))
            figure_canvas_agg.draw()
            toolbar.update()
            figure_canvas_agg.get_tk_widget().pack()

        def clearfig(figname):
            axes=figname.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_agg):
            figure_canvas_agg.get_tk_widget().forget()
            plt.close('all')

        while True:
            event6, values6 = window6.read()

            if event6 == None:
                break

            if event6 == '-confirm1-':
                try:
                    im = float(values6['-im-'])
                    dm = float(values6['-dm-'])
                    table_content.append([len(table_content)+1, im, dm])
                    window6['-table-'].update(table_content)

                except:
                    sg.popup("请输入正确的地震动强度指标IM和结构需求参数DM", title='提示', text_color='red', font=20)

            if event6 == '-delete1-':
                index = values6['-table-']
                
                if len(index) == 0:
                    sg.popup("请先在表中点击需要删除的数据行", title='提示', text_color='red', font=20)
                
                else:
                    no = -1
                    
                    for i in index:
                        no += 1
                        table_content.pop(i-no)
                    
                    for id in range(len(table_content)):
                        table_content[id][0] = id+1
                    
                    window6['-table-'].update(table_content)
            
            if event6 == '-limitstatenum-':
                choose_num = values6['-limitstatenum-']
                window6['-ls1-'].update('')
                window6['-ls2-'].update('')
                window6['-ls3-'].update('')
                window6['-ls4-'].update('')
                window6['-ls5-'].update('')
                window6['-ls6-'].update('')
                window6['-ls1query-'].update('')
                window6['-ls2query-'].update('')
                window6['-ls3query-'].update('')
                window6['-ls4query-'].update('')
                window6['-ls5query-'].update('')
                window6['-ls6query-'].update('')
                window6['-imquery-'].update('')
                window6['-equal1-'].update('ln(DM) = b + a·ln(IM)')
                window6['-nihenum-'].update('一次')
                window6['-startcoord-'].update('')
                window6['-endcoord-'].update('')
                delete_fig_agg(figure_canvas_agg1)
                delete_fig_agg(figure_canvas_agg2)
                
                if choose_num == '1':
                    window6['-textls2-'].update(visible=False)
                    window6['-ls2-'].update(visible=False)
                    window6['-textquery2-'].update(visible=False)
                    window6['-ls2query-'].update(visible=False)
                    window6['-textls3-'].update(visible=False)
                    window6['-ls3-'].update(visible=False)
                    window6['-textquery3-'].update(visible=False)
                    window6['-ls3query-'].update(visible=False)
                    window6['-textls4-'].update(visible=False)
                    window6['-ls4-'].update(visible=False)
                    window6['-textquery4-'].update(visible=False)
                    window6['-ls4query-'].update(visible=False)
                    window6['-textls5-'].update(visible=False)
                    window6['-ls5-'].update(visible=False)
                    window6['-textquery5-'].update(visible=False)
                    window6['-ls5query-'].update(visible=False)
                    window6['-textls6-'].update(visible=False)
                    window6['-ls6-'].update(visible=False)
                    window6['-textquery6-'].update(visible=False)
                    window6['-ls6query-'].update(visible=False)
                elif choose_num == '2':
                    window6['-textls2-'].update(visible=True)
                    window6['-ls2-'].update(visible=True)
                    window6['-textquery2-'].update(visible=True)
                    window6['-ls2query-'].update(visible=True)
                    window6['-textls3-'].update(visible=False)
                    window6['-ls3-'].update(visible=False)
                    window6['-textquery3-'].update(visible=False)
                    window6['-ls3query-'].update(visible=False)
                    window6['-textls4-'].update(visible=False)
                    window6['-ls4-'].update(visible=False)
                    window6['-textquery4-'].update(visible=False)
                    window6['-ls4query-'].update(visible=False)
                    window6['-textls5-'].update(visible=False)
                    window6['-ls5-'].update(visible=False)
                    window6['-textquery5-'].update(visible=False)
                    window6['-ls5query-'].update(visible=False)
                    window6['-textls6-'].update(visible=False)
                    window6['-ls6-'].update(visible=False)
                    window6['-textquery6-'].update(visible=False)
                    window6['-ls6query-'].update(visible=False)
                elif choose_num == '3':
                    window6['-textls2-'].update(visible=True)
                    window6['-ls2-'].update(visible=True)
                    window6['-textquery2-'].update(visible=True)
                    window6['-ls2query-'].update(visible=True)
                    window6['-textls3-'].update(visible=True)
                    window6['-ls3-'].update(visible=True)
                    window6['-textquery3-'].update(visible=True)
                    window6['-ls3query-'].update(visible=True)
                    window6['-textls4-'].update(visible=False)
                    window6['-ls4-'].update(visible=False)
                    window6['-textquery4-'].update(visible=False)
                    window6['-ls4query-'].update(visible=False)
                    window6['-textls5-'].update(visible=False)
                    window6['-ls5-'].update(visible=False)
                    window6['-textquery5-'].update(visible=False)
                    window6['-ls5query-'].update(visible=False)
                    window6['-textls6-'].update(visible=False)
                    window6['-ls6-'].update(visible=False)
                    window6['-textquery6-'].update(visible=False)
                    window6['-ls6query-'].update(visible=False)
                elif choose_num == '4':
                    window6['-textls2-'].update(visible=True)
                    window6['-ls2-'].update(visible=True)
                    window6['-textquery2-'].update(visible=True)
                    window6['-ls2query-'].update(visible=True)
                    window6['-textls3-'].update(visible=True)
                    window6['-ls3-'].update(visible=True)
                    window6['-textquery3-'].update(visible=True)
                    window6['-ls3query-'].update(visible=True)
                    window6['-textls4-'].update(visible=True)
                    window6['-ls4-'].update(visible=True)
                    window6['-textquery4-'].update(visible=True)
                    window6['-ls4query-'].update(visible=True)
                    window6['-textls5-'].update(visible=False)
                    window6['-ls5-'].update(visible=False)
                    window6['-textquery5-'].update(visible=False)
                    window6['-ls5query-'].update(visible=False)
                    window6['-textls6-'].update(visible=False)
                    window6['-ls6-'].update(visible=False)
                    window6['-textquery6-'].update(visible=False)
                    window6['-ls6query-'].update(visible=False)
                elif choose_num == '5':
                    window6['-textls2-'].update(visible=True)
                    window6['-ls2-'].update(visible=True)
                    window6['-textquery2-'].update(visible=True)
                    window6['-ls2query-'].update(visible=True)
                    window6['-textls3-'].update(visible=True)
                    window6['-ls3-'].update(visible=True)
                    window6['-textquery3-'].update(visible=True)
                    window6['-ls3query-'].update(visible=True)
                    window6['-textls4-'].update(visible=True)
                    window6['-ls4-'].update(visible=True)
                    window6['-textquery4-'].update(visible=True)
                    window6['-ls4query-'].update(visible=True)
                    window6['-textls5-'].update(visible=True)
                    window6['-ls5-'].update(visible=True)
                    window6['-textquery5-'].update(visible=True)
                    window6['-ls5query-'].update(visible=True)
                    window6['-textls6-'].update(visible=False)
                    window6['-ls6-'].update(visible=False)
                    window6['-textquery6-'].update(visible=False)
                    window6['-ls6query-'].update(visible=False)
                elif choose_num == '6':
                    window6['-textls2-'].update(visible=True)
                    window6['-ls2-'].update(visible=True)
                    window6['-textquery2-'].update(visible=True)
                    window6['-ls2query-'].update(visible=True)
                    window6['-textls3-'].update(visible=True)
                    window6['-ls3-'].update(visible=True)
                    window6['-textquery3-'].update(visible=True)
                    window6['-ls3query-'].update(visible=True)
                    window6['-textls4-'].update(visible=True)
                    window6['-ls4-'].update(visible=True)
                    window6['-textquery4-'].update(visible=True)
                    window6['-ls4query-'].update(visible=True)
                    window6['-textls5-'].update(visible=True)
                    window6['-ls5-'].update(visible=True)
                    window6['-textquery5-'].update(visible=True)
                    window6['-ls5query-'].update(visible=True)
                    window6['-textls6-'].update(visible=True)
                    window6['-ls6-'].update(visible=True)
                    window6['-textquery6-'].update(visible=True)
                    window6['-ls6query-'].update(visible=True)

            if event6 == '-nihenum-':
                if values6['-nihenum-'] == "二次":
                    window6['-equal1-'].update('ln(DM) = c + b·ln(IM) + a·ln(IM)^2')
                else:
                    window6['-equal1-'].update('ln(DM) = b + a·ln(IM)')

            if event6 == '-confirm2-':
                try:
                    clearfig(fig1)
                    clearfig(fig2)
                    DIs=[]
                    choose_num = values6['-limitstatenum-']
                    ims_start = float(values6['-startcoord-'])
                    ims_end = float(values6['-endcoord-'])
                    ims_fragilitycurve = np.linspace(ims_start, ims_end, 500)
                    IMs, DMs = [], []
                    for i in range(len(table_content)):
                        IMs.append(table_content[i][1])
                        DMs.append(table_content[i][2])
                    ln_IMs=np.log(IMs)
                    ln_DMs=np.log(DMs)
                    scatterfig(fig1, figure_canvas_agg1, toolbar1, ln_IMs, ln_DMs)

                    if values6['-nihenum-'] == '一次':
                        deg = 1
                        a, b = fragility_ida(ln_IMs, ln_DMs, deg)
                        window6['-equal1-'].update(f"ln(DM)={round(b, 3)}+{round(a, 3)}*ln(IM)")
                        im_min, im_max = np.min(ln_IMs), np.max(ln_IMs)
                        ims_nihe = np.linspace(im_min, im_max, 500)
                        dms_nihe = a * ims_nihe + b
                        updatefig(fig1, figure_canvas_agg1, toolbar1, ims_nihe, dms_nihe, 1, '拟合公式')

                        if choose_num == '1':
                            DI1 = float(values6['-ls1-'])
                            DIs.append(DI1)
                            for i in range(len(DIs)):
                                zs = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')
                        
                        if choose_num == '2':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            for i in range(len(DIs)):
                                zs = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '3':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            for i in range(len(DIs)):
                                zs = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '4':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DI4 = float(values6['-ls4-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            DIs.append(DI4)
                            for i in range(len(DIs)):
                                zs = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '5':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DI4 = float(values6['-ls4-'])
                            DI5 = float(values6['-ls5-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            DIs.append(DI4)
                            DIs.append(DI5)
                            for i in range(len(DIs)):
                                zs = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '6':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DI4 = float(values6['-ls4-'])
                            DI5 = float(values6['-ls5-'])
                            DI6 = float(values6['-ls6-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            DIs.append(DI4)
                            DIs.append(DI5)
                            DIs.append(DI6)
                            for i in range(len(DIs)):
                                zs = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')
                        

                    if values6['-nihenum-'] == '二次':
                        deg = 2
                        a, b, c = fragility_ida(ln_IMs, ln_DMs, deg)
                        window6['-equal1-'].update(f"ln(DM)={round(c, 3)}+{round(b, 3)}*ln(IM)+{round(a, 3)}*ln(IM)^2")
                        im_min, im_max = np.min(ln_IMs), np.max(ln_IMs)
                        ims_nihe = np.linspace(im_min, im_max, 500)
                        ims_nihe_erci = []
                        for im_nihe_erci in ims_nihe:
                            ims_nihe_erci.append(pow(im_nihe_erci, 2))
                        ims_nihe_erci = np.array(ims_nihe_erci)
                        dms_nihe = a * ims_nihe_erci + b * ims_nihe + c
                        updatefig(fig1, figure_canvas_agg1, toolbar1, ims_nihe, dms_nihe, 1, '拟合公式')


                        if choose_num == '1':
                            DI1 = float(values6['-ls1-'])
                            DIs.append(DI1)
                            ims_fragilitycurve_erci = []
                            for im_erci in ims_fragilitycurve:
                                ims_fragilitycurve_erci.append(pow(np.log(im_erci), 2))
                            ims_fragilitycurve_erci = np.array(ims_fragilitycurve_erci)
                            for i in range(len(DIs)):
                                zs = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')
                        
                        if choose_num == '2':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            ims_fragilitycurve_erci = []
                            for im_erci in ims_fragilitycurve:
                                ims_fragilitycurve_erci.append(pow(np.log(im_erci), 2))
                            ims_fragilitycurve_erci = np.array(ims_fragilitycurve_erci)
                            for i in range(len(DIs)):
                                zs = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '3':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            ims_fragilitycurve_erci = []
                            for im_erci in ims_fragilitycurve:
                                ims_fragilitycurve_erci.append(pow(np.log(im_erci), 2))
                            ims_fragilitycurve_erci = np.array(ims_fragilitycurve_erci)
                            for i in range(len(DIs)):
                                zs = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '4':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DI4 = float(values6['-ls4-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            DIs.append(DI4)
                            ims_fragilitycurve_erci = []
                            for im_erci in ims_fragilitycurve:
                                ims_fragilitycurve_erci.append(pow(np.log(im_erci), 2))
                            ims_fragilitycurve_erci = np.array(ims_fragilitycurve_erci)
                            for i in range(len(DIs)):
                                zs = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '5':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DI4 = float(values6['-ls4-'])
                            DI5 = float(values6['-ls5-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            DIs.append(DI4)
                            DIs.append(DI5)
                            ims_fragilitycurve_erci = []
                            for im_erci in ims_fragilitycurve:
                                ims_fragilitycurve_erci.append(pow(np.log(im_erci), 2))
                            ims_fragilitycurve_erci = np.array(ims_fragilitycurve_erci)
                            for i in range(len(DIs)):
                                zs = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                        if choose_num == '6':
                            DI1 = float(values6['-ls1-'])
                            DI2 = float(values6['-ls2-'])
                            DI3 = float(values6['-ls3-'])
                            DI4 = float(values6['-ls4-'])
                            DI5 = float(values6['-ls5-'])
                            DI6 = float(values6['-ls6-'])
                            DIs.append(DI1)
                            DIs.append(DI2)
                            DIs.append(DI3)
                            DIs.append(DI4)
                            DIs.append(DI5)
                            DIs.append(DI6)
                            ims_fragilitycurve_erci = []
                            for im_erci in ims_fragilitycurve:
                                ims_fragilitycurve_erci.append(pow(np.log(im_erci), 2))
                            ims_fragilitycurve_erci = np.array(ims_fragilitycurve_erci)
                            for i in range(len(DIs)):
                                zs = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                                Pfs = st.norm.cdf(zs)*100
                                drawfig(fig2, figure_canvas_agg2, toolbar2, ims_fragilitycurve, Pfs, 1, label=f'LS{i+1}')

                except:
                    sg.popup('请检查输入信息是否完整,尤其是易损性曲线横坐标的设置！', title='提示', text_color='red', font=20)

            if event6 == '-query-':
                try:
                    aim_im = float(values6['-imquery-'])
                    
                    if values6['-nihenum-'] == '一次':
                        for i in range(len(DIs)):
                            zs_i = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                            Pfs_i = st.norm.cdf(zs_i)*100
                            for id in range(len(ims_fragilitycurve)):
                                if ims_fragilitycurve[id]< aim_im <=ims_fragilitycurve[id+1]:
                                    Pf_query = (Pfs_i[id+1]-Pfs_i[id])/(ims_fragilitycurve[id+1]-ims_fragilitycurve[id])*(aim_im-ims_fragilitycurve[id])+Pfs_i[id]
                            window6[f'-ls{i+1}query-'].update(f'{round(Pf_query, 2)}%')
                        
                    if values6['-nihenum-'] == '二次':
                        for i in range(len(DIs)):
                            zs_i = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                            Pfs_i = st.norm.cdf(zs_i)*100
                            for id in range(len(ims_fragilitycurve)):
                                if ims_fragilitycurve[id]< aim_im <=ims_fragilitycurve[id+1]:
                                    Pf_query = (Pfs_i[id+1]-Pfs_i[id])/(ims_fragilitycurve[id+1]-ims_fragilitycurve[id])*(aim_im-ims_fragilitycurve[id])+Pfs_i[id]
                            window6[f'-ls{i+1}query-'].update(f'{round(Pf_query, 2)}%')

                except:
                    sg.popup('请检查输入信息是否完整', title='提示', text_color='red', font=20)

            if event6 == '-save-':
                savefile_location = values6['-savelocation-']
                
                try:
                    with open(savefile_location, 'w') as f:
                        
                        if values6['-nihenum-'] == '一次':
                            for i in range(len(DIs)):
                                f.writelines(f'损伤状态 LS{i+1}易损性曲线:\n')
                                f.writelines('IM Pf\n')
                                zs_i = (b + a * np.log(ims_fragilitycurve) - np.log(DIs[i]))/0.5
                                Pfs_i = st.norm.cdf(zs_i)*100
                                for im, p in zip(ims_fragilitycurve, Pfs_i):
                                    f.writelines('%f %f\n' %(im, p))
                                f.writelines('\n')
                        
                        if values6['-nihenum-'] == '二次':
                            for i in range(len(DIs)):
                                f.writelines(f'损伤状态 LS{i+1}易损性曲线:\n')
                                f.writelines('IM Pf\n')
                                zs_i = (c + b * np.log(ims_fragilitycurve) + a * ims_fragilitycurve_erci - np.log(DIs[i]))/0.5
                                Pfs_i = st.norm.cdf(zs_i)*100
                                for im, p in zip(ims_fragilitycurve, Pfs_i):
                                    f.writelines('%f %f\n' %(im, p))
                                f.writelines('\n')
                    
                    sg.popup('保存成功', title='提示', text_color='green', font=20)
                
                except:
                    sg.popup('请先选择保存文件名', title='提示', text_color='red', font=20)

            if event6 == '-clearall-':
                window6['-im-'].update('')
                window6['-dm-'].update('')
                table_content = []
                window6['-table-'].update('')
                window6['-ls1-'].update('')
                window6['-ls2-'].update('')
                window6['-ls3-'].update('')
                window6['-ls4-'].update('')
                window6['-ls5-'].update('')
                window6['-ls6-'].update('')
                window6['-ls1query-'].update('')
                window6['-ls2query-'].update('')
                window6['-ls3query-'].update('')
                window6['-ls4query-'].update('')
                window6['-ls5query-'].update('')
                window6['-ls6query-'].update('')
                window6['-imquery-'].update('')
                window6['-equal1-'].update('ln(DM) = b + a·ln(IM)')
                window6['-nihenum-'].update('一次')
                window6['-limitstatenum-'].update('4')
                delete_fig_agg(figure_canvas_agg1)
                delete_fig_agg(figure_canvas_agg2)
                window6['-startcoord-'].update('')
                window6['-endcoord-'].update('')
                window6['-savelocation-'].update('')
                window6['-textls2-'].update(visible=True)
                window6['-ls2-'].update(visible=True)
                window6['-textquery2-'].update(visible=True)
                window6['-ls2query-'].update(visible=True)
                window6['-textls3-'].update(visible=True)
                window6['-ls3-'].update(visible=True)
                window6['-textquery3-'].update(visible=True)
                window6['-ls3query-'].update(visible=True)
                window6['-textls4-'].update(visible=True)
                window6['-ls4-'].update(visible=True)
                window6['-textquery4-'].update(visible=True)
                window6['-ls4query-'].update(visible=True)
                window6['-textls5-'].update(visible=False)
                window6['-ls5-'].update(visible=False)
                window6['-textquery5-'].update(visible=False)
                window6['-ls5query-'].update(visible=False)
                window6['-textls6-'].update(visible=False)
                window6['-ls6-'].update(visible=False)
                window6['-textquery6-'].update(visible=False)
                window6['-ls6query-'].update(visible=False)
                sg.popup('清除成功！', font=20)

            if event6 == sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
                    window6.close()
                    window0.close()
    
            if event6 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提醒', text_color='red', font=20)=="Yes":
                    window6.close()
                    win6_active=False
                    window0.UnHide()
                


# ## -----------window7----------window7----------window7----------window7----------window7----------window7----------window7----------  
    if event0 == "-skeletoncurve-" and not win7_active:
        win7_active = True
        window0.Hide()

        table_content = []
        total_filenames_location = []

        layout71 = [
            [sg.FilesBrowse("添加滞回曲线文件", key='-openfiles-', font=20, file_types=((('ALL Files', '*.txt'),)), target='-filenames-'), sg.In("", key='-filenames-', font=20), sg.B('添加', key='-confirm1-', font=20), sg.B("清除选择", key="-clearchoose-", font=20), sg.B("删除", key="-delete-", font=20), sg.Text('数据划分点数', font=20), sg.In('500', key='-pointnum-', font=20, size=(6, 1)), sg.B('绘图', key="-draw-", font=20)],
            [sg.Table(headings=['序号', '文件名称'],
                    values=table_content,
                    expand_x=True,
                    expand_y=True, 
                    font=20,
                    key='-table-',
                    justification='center',
                    select_mode="extended",
                    enable_events=True)],
            [sg.Canvas(key='-toolbar-', size=(80, 1))],
            [sg.Canvas(key='-canvas-', expand_x=True, expand_y=True)]
        ]

        layout72 = [
            [sg.FolderBrowse('选择保存文件夹', target='-savelocation-', font=20), sg.In('', key='-savelocation-', visible=True, font=20), sg.B('保存', key='-save-', font=20), sg.B('清除所有', key='-clearall-', font=20), sg.B('返回至初始界面', key='-back-', button_color=('red', 'skyblue'), font=20)]
        ]

        layout70 = [
            [sg.Frame('结构滞回曲线选择', layout71, title_color='red', font=('楷体'), expand_x=True, element_justification='center', expand_y=True)],
            [sg.Frame('批量导出骨架曲线文件', layout72, title_color='red', font=('楷体'), expand_x=True, element_justification='center', expand_y=True)]
        ]

        layout7 = [[sg.Column(layout70, scrollable=True, vertical_scroll_only=True, size=(1400, 1050))]]
        window7 = sg.Window('土木工程数据常用处理软件 ----- 骨架曲线提取', layout7, location=(300, 10), finalize=True, resizable=True, enable_close_attempted_event=True)

        def exractenvelop(disps, forces, point_num):
            disp_max = np.max(disps)
            xs = np.linspace(0, disp_max, point_num)
            ys = []

            for i in range(len(xs)-1):
                yi = []
                for id in range(len(disps)-1):
                    if disps[id]<=xs[i]<disps[id+1]:
                        zanshi_y = (forces[id+1]-forces[id])/(disps[id+1]-disps[id])*(xs[i]-disps[id])+forces[id]
                        yi.append(zanshi_y)
                yi = np.array(yi)
                yi_max = np.max(yi)
                ys.append(yi_max)
            xs_envelope = xs[ :-1]
            ys_envelope = ys

            for id in range(len(disps)):
                for i in range(len(xs_envelope)-1):
                    if xs_envelope[i]<=disps[id]<xs_envelope[i+1]:
                        zanshi_y = (ys_envelope[i+1]-ys_envelope[i])/(xs_envelope[i+1]-xs_envelope[i])*(disps[id]-xs_envelope[i])+ys_envelope[i]
                        delta = zanshi_y - forces[id]
                if delta <= 0:
                    id_chaoyue = id
                    break
                else:
                    continue

            for i in range(1, len(xs_envelope)):
                if xs_envelope[i-1]<disps[id_chaoyue]<=xs_envelope[i]:
                    i_chaoyue = i

            new_disps_after = np.array(xs_envelope[i_chaoyue:])
            new_forces_after = np.array(ys_envelope[i_chaoyue:])

            no = 0
            while True:
                no += 1
                if forces[id_chaoyue-no]>=0 and forces[id_chaoyue-no-1]<0:
                    id_start = id_chaoyue-no
                    break
                if id_chaoyue-no == 0:
                    id_start = 0
                    break

            new_disps_before = np.array(disps[id_start:id_chaoyue])
            new_forces_before = np.array(forces[id_start:id_chaoyue])

            new_disps_ = np.concatenate((new_disps_before, new_disps_after), axis=0)
            new_forces_ = np.concatenate((new_forces_before, new_forces_after), axis=0)

            return new_disps_, new_forces_

        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)
        
        fig = matplotlib.figure.Figure(figsize=(8, 5.7))
        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False
        fig.add_subplot(111)
        axes = fig.axes
        axes[0].set_xlabel('位移')
        axes[0].set_ylabel('荷载')
        figure_canvas_agg = FigureCanvasTkAgg(fig, window7['-canvas-'].TKCanvas)
        toolbar = Toolbar(figure_canvas_agg, window7['-toolbar-'].TKCanvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack()

        def updatefig(disps, forces, alpha, label):
            axes=fig.axes
            axes[0].set_xlabel('位移')
            axes[0].set_ylabel('荷载')
            axes[0].plot(disps, forces, alpha=alpha, label=label)
            axes[0].legend()
            toolbar.update()
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def scatterfig(point_disp, point_force):
            axes=fig.axes
            axes[0].scatter(point_disp, point_force, c='black', zorder=2)
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def clearfig():
            axes=fig.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_agg):
            figure_canvas_agg.get_tk_widget().forget()
            plt.close('all')

        while True:
            event7, values7 = window7.read()

            if event7 == None:
                break

            if event7 == '-confirm1-':
                try:
                    window7['-filenames-'].update('')
                    filenames = str(values7['-filenames-']).split(';')
                    current_filenames = []
                    for i in range(len(table_content)):
                        current_filenames.append(table_content[i][1])
                    for filename in filenames:
                        dataname = filename.split('/')[-1]
                        if dataname in current_filenames:
                            sg.popup('请勿重复添加文件', title='提示', text_color='red')
                            break
                        else:
                            total_filenames_location.append(filename)
                            table_content.append([len(table_content)+1, dataname])
                            window7['-table-'].update(table_content)
                except:
                    sg.popup('请先选择要打开的文件', title='提示', text_color='red', font=20)

            if event7 == '-clearchoose-':
                window7["-filenames-"].update("")

            if event7 == '-delete-':
                index = values7["-table-"]
                num = -1
                for i in index:
                    num += 1
                    id = i - num
                    table_content.pop(id)
                    total_filenames_location.pop(id)
                for id in range(len(table_content)):
                    table_content[id][0] = id+1
                window7['-table-'].update(table_content)

            if event7 == "-draw-":
                try:
                    clearfig()
                    pointnum = int(values7['-pointnum-'])
                    index = values7['-table-']
                    
                    if len(index)==0:
                        sg.popup('请先在表中点击一组需要绘图的数据', title='提示', text_color='red', font=20)
                    elif len(index)>1:
                        sg.popup('每次仅支持查看一组数据', title='提示', text_color='red', font=20)
                    else:
                        choose_filename_location = total_filenames_location[index[0]]
                        datas = np.loadtxt(choose_filename_location, dtype=float, unpack=True)
                        disps = datas[0]
                        forces = datas[1]
                        new_disps_, new_forces_ = exractenvelop(disps, forces, pointnum)
                        updatefig(disps, forces, 0.8, label='滞回曲线')
                        updatefig(new_disps_, new_forces_, 1, label='骨架曲线')
                except:
                    sg.popup('请检查是否完整填写信息', title='提示', text_color='red', font=20)

            if event7 == "-save-":
                try:
                    folder_location = values7['-savelocation-']
                    pointnum = int(values7['-pointnum-'])

                    for i in range(len(table_content)):
                        filename = table_content[i][1]
                        file_location = total_filenames_location[i]
                        savefile_location = folder_location+'/骨架曲线_'+filename
                        datas = np.loadtxt(file_location, dtype=float, unpack=True)
                        disps = datas[0]
                        forces = datas[1]
                        new_disps_, new_forces_ = exractenvelop(disps, forces, pointnum)

                        with open(savefile_location, 'w') as f:
                            for disp, force in zip(new_disps_, new_forces_):
                                f.writelines("%f %f\n" %(disp, force))
                    sg.popup('保存成功', title='提示', text_color='green', font=20)
                except:
                    sg.popup('请先选择保存文件夹', title='提示', text_color='red', font=20)

            if event7 == "-clearall-":
                window7['-filenames-'].update('')
                table_content = []
                total_filenames_location = []
                window7['-table-'].update(table_content)
                window7['-pointnum-'].update('500')
                delete_fig_agg(figure_canvas_agg)
                window7['-savelocation-'].update('')

            if event7==sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
                    window7.close()
                    window0.close()

            if event7 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提醒', text_color='red', font=20)=="Yes":
                    window7.close()
                    win7_active=False
                    window0.UnHide()

            

## -----------window8----------window8----------window8----------window8----------window8----------window8----------window8----------
    if event0 == "-curvefilter-" and not win8_active:
        win8_active = True
        window0.Hide()

        table_content = []
        total_filenames_location = []

        layout81 = [
            [sg.FilesBrowse("添加试验曲线文件", key='-openfiles-', font=20, file_types=((('ALL Files', '*.txt'),)), target='-filenames-'), sg.In("", key='-filenames-', font=20), sg.B('添加', key='-confirm1-', font=20), sg.B("清除选择", key="-clearchoose-", font=20), sg.B("删除", key="-delete-", font=20), sg.B('选择并绘图', key="-draw-", font=20)],
            [sg.Text('卷积核长度', font=20), sg.Slider((1, 201), default_value=1, tick_interval=20, border_width=5, key='-jjlen-', font=20, enable_events=True, orientation='horizontal', expand_x=True)],
            [sg.Table(headings=['序号', '文件名称', '使用卷积核长度'],
                    values=table_content,
                    expand_x=True,
                    expand_y=True, 
                    font=20,
                    key='-table-',
                    num_rows=6,
                    justification='center',
                    select_mode="extended",
                    enable_events=True)],
            [sg.Canvas(key='-toolbar-', size=(80, 1))],
            [sg.Canvas(key='-canvas-', expand_x=True, expand_y=True)]
        ]

        layout82 = [
            [sg.FolderBrowse('选择保存文件夹', target='-savelocation-', font=20), sg.In('', key='-savelocation-', visible=True, font=20), sg.B('保存', key='-save-', font=20), sg.B('清除所有', key='-clearall-', font=20), sg.B('返回至初始界面', key='-back-', button_color=('red', 'skyblue'), font=20)]
        ]

        layout80 = [
            [sg.Frame('试验曲线选择', layout81, title_color='red', font=('楷体'), expand_x=True, element_justification='center', expand_y=True)],
            [sg.Frame('批量导出文件', layout82, title_color='red', font=('楷体'), expand_x=True, element_justification='center', expand_y=True)]
        ]

        layout8 = [[sg.Column(layout80, scrollable=True, vertical_scroll_only=True, size=(1200, 1100))]]

        window8 = sg.Window('土木工程数据常用处理软件 ----- 试验曲线降噪', layout8, location=(300, 10), finalize=True, resizable=True, enable_close_attempted_event=True)
        
        class Toolbar(NavigationToolbar2Tk):
            def __init__(self, *args, **kwargs):
                super(Toolbar, self).__init__(*args, **kwargs)

        def moving_average(interval, windowsize):
            window = np.ones(int(windowsize)) / float(windowsize)
            re = np.convolve(interval, window, 'same')
            return re

        fig = matplotlib.figure.Figure(figsize=(8, 5.7))
        plt.rcParams['font.sans-serif'] = ['Microsoft Yahei']
        plt.rcParams['axes.unicode_minus'] = False
        fig.add_subplot(111)
        axes = fig.axes
        figure_canvas_agg = FigureCanvasTkAgg(fig, window8['-canvas-'].TKCanvas)
        toolbar = Toolbar(figure_canvas_agg, window8['-toolbar-'].TKCanvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack()

        def updatefig(disps, forces, alpha, label):
            axes=fig.axes
            axes[0].plot(disps, forces, alpha=alpha, label=label)
            axes[0].legend()
            toolbar.update()
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def scatterfig(point_disp, point_force):
            axes=fig.axes
            axes[0].scatter(point_disp, point_force, c='black', zorder=2)
            figure_canvas_agg.draw()
            figure_canvas_agg.get_tk_widget().pack()

        def clearfig():
            axes=fig.axes
            axes[0].cla()

        def delete_fig_agg(figure_canvas_agg):
            figure_canvas_agg.get_tk_widget().forget()
            plt.close('all')


        while True:
            event8, values8 = window8.read()

            if event8 == None:
                break

            if event8 == '-confirm1-':
                try:
                    window8['-filenames-'].update('')
                    filenames = str(values8['-filenames-']).split(';')
                    current_filenames = []
                    for i in range(len(table_content)):
                        current_filenames.append(table_content[i][1])
                    for filename in filenames:
                        dataname = filename.split('/')[-1]
                        if dataname in current_filenames:
                            sg.popup('请勿重复添加文件', title='提示', text_color='red')
                            break
                        else:
                            total_filenames_location.append(filename)
                            table_content.append([len(table_content)+1, dataname, 1])
                            window8['-table-'].update(table_content)
                except:
                    sg.popup('请先选择要打开的文件', title='提示', text_color='red', font=20)

            if event8 == '-clearchoose-':
                window8["-filenames-"].update("")

            if event8 == '-delete-':
                index = values8["-table-"]
                num = -1
                for i in index:
                    num += 1
                    id = i - num
                    table_content.pop(id)
                    total_filenames_location.pop(id)
                for id in range(len(table_content)):
                    table_content[id][0] = id+1
                window8['-table-'].update(table_content)

            if event8 == '-draw-':
                clearfig()
                index = values8['-table-']
                try:
                    if len(index) ==0:
                        sg.popup('请先在表格中选择一组数据', title='提示', text_color='red', font=20)
                    elif len(index) >1:
                        sg.popup('每次仅支持查看一组数据', title='提示', text_color='red', font=20)
                    else:
                        filename_location = total_filenames_location[index[0]]
                        no = table_content[index[0]][0]
                        datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                        disps = datas[0]
                        forces = datas[1]
                        updatefig(disps, forces, 0.8, label=f'试验曲线{no}')
                except:
                    sg.popup('操作失误', title='提示', text_color='red', font=20)

            if event8 == '-jjlen-':
                try:
                    clearfig()
                    juanji_len = int(values8['-jjlen-'])
                    filename_location = total_filenames_location[index[0]]
                    datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                    disps = datas[0]
                    forces = datas[1]
                    windowsize = table_content[index[0]][2] = juanji_len
                    no = table_content[index[0]][0]
                    window8['-table-'].update(table_content)
                    updatefig(disps, forces, 0.8, label=f'试验曲线{no}')
                    forces_filter = moving_average(forces, windowsize)
                    updatefig(disps, forces_filter, 1, label='降噪后曲线')  
                except:
                    sg.popup('请先在表格中选择一组数据并点击确定绘图', title='提示', text_color='red', font=20)
                    window8['-jjlen-'].update('1')
            
            if event8 == '-save-':
                try:
                    for i in range(len(table_content)):
                        folder_location = values8['-savelocation-']
                        filename_location = total_filenames_location[i]
                        savefilename = table_content[i][1]
                        savefile_location = folder_location+'/降噪后_'+savefilename
                        windowsize = int(table_content[i][2])
                        datas = np.loadtxt(filename_location, dtype=float, unpack=True)
                        disps = datas[0]
                        forces = datas[1]
                        forces_filter = moving_average(forces, windowsize)

                        with open(savefile_location, 'w') as f:
                            for d, force in zip(disps, forces_filter):
                                f.writelines("%f %f\n" %(d, force))
                    sg.popup('保存成功', title='提示', text_color='green', font=20)
                
                except:
                    sg.popup('请先选择保存文件夹', title='提示', text_color='red', font=20)
            
            if event8 == '-clearall-':
                window8['-filenames-'].update('')
                table_content = []
                total_filenames_location = []
                window8['-table-'].update(table_content)
                window8['-jjlen-'].update('1')
                window8['-savelocation-'].update('')
                delete_fig_agg(figure_canvas_agg)
                sg.popup('清除成功', title='提示', text_color='green', font=20)

            if event8==sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
                if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
                    window8.close()
                    window0.close()

            if event8 == '-back-':
                if sg.popup_yes_no("请问您是否要回到初始界面", title='提醒', text_color='red', font=20)=="Yes":
                    window8.close()
                    win8_active=False
                    window0.UnHide()

            











                    



                    

                    
    if event0 == sg.WINDOW_CLOSE_ATTEMPTED_EVENT:
        if sg.popup_yes_no("请问您是否要退出程序", title='提醒', text_color='red', font=20)=="Yes":
            window0.close()



