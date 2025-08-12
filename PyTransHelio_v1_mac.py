"""
v6.0 修复了选区跨越0度经线时的计算错误


"""
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename
import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
from matplotlib.widgets import PolygonSelector
from sunpy.map import Map
from sunpy.coordinates import HeliographicStonyhurst
from sunpy.coordinates import HeliographicCarrington
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits

#添加内容，强制程序使用TkAgg终端
import matplotlib
matplotlib.use('TkAgg')




class FITSViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("PyTransHelio")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)  # 绑定关闭事件

        # 数据属性初始化
        self.data = None
        self.im = None
        self.colorbar = None
        self.rect = None
        self.hdul = None
        self.smap = None
        self.current_filename = None
        self.txt_path_var = tk.StringVar()

        # 配置三列布局
        root.grid_columnconfigure(0, weight=2)  # 图像列
        root.grid_columnconfigure(1, weight=1)  # 控制列
        root.grid_columnconfigure(2, weight=2)  # 输出列

        # === 图像列 ===
        img_frame = tk.Frame(root)
        img_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

        # 调整图形布局为垂直排列
        self.fig = plt.figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111)  # 主图区域
        self.canvas = FigureCanvasTkAgg(self.fig, master=img_frame)

        # matplotlib工具栏
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
        self.toolbar = NavigationToolbar2Tk(self.canvas, img_frame)
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # === 控制列 ===
        ctrl_frame = tk.Frame(root)
        ctrl_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)

        # 创建控制列组件
        self._create_file_controls(ctrl_frame)  # 行1：文件操作按钮
        self._create_vmin_vmax_controls(ctrl_frame)  # 第二行：vmin/vmax
        self._create_group_selector(ctrl_frame)  # 第三行：组选择
        self._create_n_settings(ctrl_frame)  # 第四行：N足点设置
        self._create_s_settings(ctrl_frame)  # 第五行：S足点设置
        self._create_calculate_button(ctrl_frame)  # 第六行：计算按钮
        self._create_txt_controls(ctrl_frame)  # 第七-九行：TXT相关控件

        # === 输出列 ===
        out_frame = tk.Frame(root)
        out_frame.grid(row=0, column=2, sticky="nsew", padx=5, pady=5)

        # 当前文件名显示
        self.filename_label = tk.Label(out_frame, text='current file：not selected', justify=tk.LEFT, anchor='w')
        self.filename_label.pack(pady=2, fill=tk.X)

        self.group_data = {
            "group1": None,
            "group2": None
        }

        # 观测时间显示
        self.obs_time_label = tk.Label(out_frame, text="observation time：not selected", justify=tk.LEFT, anchor='w')
        self.obs_time_label.pack(pady=2, fill=tk.X)

        # 点击坐标显示
        self.coordinate_label = tk.Label(out_frame, text="coordination：none", justify=tk.LEFT, anchor='w')
        self.coordinate_label.pack(pady=2, fill=tk.X)

        # 中心坐标点显示
        self.center_coord_label = tk.Label(out_frame, text="carrington：none", justify=tk.LEFT, anchor='w')
        self.center_coord_label.pack(pady=2, fill=tk.X)

        # 输出信息显示（如报错信息等）
        self.output_label = tk.Label(out_frame, text="", justify=tk.LEFT, fg="red", anchor='w')
        self.output_label.pack(pady=2, fill=tk.X)

        # 分组信息显示框架
        result_frame = tk.Frame(out_frame)
        result_frame.pack(pady=5, fill=tk.BOTH, expand=True)

        # 第一组（N足点）信息显示
        self.group1_info = tk.Label(result_frame, text="none", justify=tk.LEFT, wraplength=200)
        self.group1_info.pack(pady=2, fill=tk.X, anchor='w')

        # 第二组（S足点）信息显示
        self.group2_info = tk.Label(result_frame, text="none", justify=tk.LEFT, wraplength=200)
        self.group2_info.pack(pady=2, fill=tk.X, anchor='w')

        # 距离信息显示
        self.distance_label = tk.Label(result_frame, text="none", justify=tk.LEFT, wraplength=200)
        self.distance_label.pack(pady=2, fill=tk.X, anchor='w')

        # 初始化阈值
        self._update_polarity_group("n")
        self._update_polarity_group("s")

        # 初始化可见性
        self._update_settings_visibility()

    # ==================控件生成========================================
    def _create_file_controls(self, parent):
        """创建文件操作按钮（第一行）"""
        btn_frame = tk.Frame(parent)
        btn_frame.grid(row=0, column=0, pady=5, sticky="ew")

        self.select_button = tk.Button(btn_frame, text="Select the FITS file", command=self.select_file)
        self.select_button.pack(side=tk.LEFT, padx=2)

        self.close_button = tk.Button(btn_frame, text="Close Image", command=self._clear_image)
        self.close_button.pack(side=tk.LEFT, padx=2)

    def _create_vmin_vmax_controls(self, parent):
        """创建显示范围设置控件（第二行）"""
        v_frame = tk.Frame(parent)
        v_frame.grid(row=1, column=0, pady=5, sticky="ew")

        # vmin 控件
        self.vmin_var = tk.StringVar(value="-200")
        self.vmin_spinbox = ttk.Spinbox(v_frame, from_=-400, to=400, increment=33,
                                        width=8, textvariable=self.vmin_var)
        self.vmin_var.trace_add("write", lambda *args: self.update_image())

        # vmax 控件
        self.vmax_var = tk.StringVar(value="200")
        self.vmax_spinbox = ttk.Spinbox(v_frame, from_=-400, to=400, increment=33,
                                        width=8, textvariable=self.vmax_var)
        self.vmax_var.trace_add("write", lambda *args: self.update_image())

        # 布局
        self.vmin_spinbox.pack(side=tk.LEFT, padx=2)
        self.vmax_spinbox.pack(side=tk.LEFT, padx=2)

    def _create_group_selector(self, parent):
        """创建组选择单选框（第三行）"""
        self.group_var = tk.StringVar(value="group1")
        frame = tk.Frame(parent)
        frame.grid(row=3, column=0, pady=5, sticky="ew")

        ttk.Radiobutton(frame, text="N FootPoint", variable=self.group_var,
                        value="group1", command=self._update_settings_visibility).pack(side=tk.LEFT)
        ttk.Radiobutton(frame, text="S FootPoint", variable=self.group_var,
                        value="group2", command=self._update_settings_visibility).pack(side=tk.LEFT)
        ttk.Radiobutton(frame, text="unselected", variable=self.group_var,
                        value="none", command=self._update_settings_visibility).pack(side=tk.LEFT)

    def _create_n_settings(self, parent):
        """N足点设置（第四行）"""
        self.n_frame = tk.Frame(parent)
        self.n_frame.grid(row=4, column=0, pady=5, sticky="ew")

        # 隐藏的阈值输入框
        self.n_lower_threshold_entry = tk.Entry(self.n_frame)
        self.n_upper_threshold_entry = tk.Entry(self.n_frame)
        self.n_lower_threshold_entry.pack_forget()
        self.n_upper_threshold_entry.pack_forget()

        # 添加组标签
        tk.Label(self.n_frame, text="N Footpoint Setting").pack(side=tk.LEFT)  # 新增标签

        # 极性选择
        self.n_polarity_var = tk.StringVar(value="positive")
        ttk.Radiobutton(self.n_frame, text="positive", variable=self.n_polarity_var,
                        value="positive", command=lambda: self._update_polarity_group("n")).pack(side=tk.LEFT)
        ttk.Radiobutton(self.n_frame, text="negative", variable=self.n_polarity_var,
                        value="negative", command=lambda: self._update_polarity_group("n")).pack(side=tk.LEFT)

        # 活动区编号
        tk.Label(self.n_frame, text="AR:").pack(side=tk.LEFT)  # 新增标签
        self.n_ar_entry = ttk.Entry(self.n_frame, width=8)
        self.n_ar_entry.pack(side=tk.LEFT, padx=2)
        self.n_ar_entry.insert(0, "0")

    def _create_s_settings(self, parent):
        """S足点设置（第五行）"""
        self.s_frame = tk.Frame(parent)
        self.s_frame.grid(row=5, column=0, pady=5, sticky="ew")

        # 阈值输入框
        self.s_lower_threshold_entry = tk.Entry(self.s_frame)
        self.s_upper_threshold_entry = tk.Entry(self.s_frame)
        self.s_lower_threshold_entry.pack_forget()
        self.s_upper_threshold_entry.pack_forget()

        # 添加组标签
        tk.Label(self.s_frame, text="N Footpoint Setting").pack(side=tk.LEFT)  # 新增标签

        # 极性选择
        self.s_polarity_var = tk.StringVar(value="positive")
        ttk.Radiobutton(self.s_frame, text="positive", variable=self.s_polarity_var,
                        value="positive", command=lambda: self._update_polarity_group("s")).pack(side=tk.LEFT)
        ttk.Radiobutton(self.s_frame, text="negative", variable=self.s_polarity_var,
                        value="negative", command=lambda: self._update_polarity_group("s")).pack(side=tk.LEFT)

        # 活动区编号
        tk.Label(self.s_frame, text="AR:").pack(side=tk.LEFT)  # 新增标签
        self.s_ar_entry = ttk.Entry(self.s_frame, width=8)
        self.s_ar_entry.pack(side=tk.LEFT, padx=2)
        self.s_ar_entry.insert(0, "0")

    def _create_calculate_button(self, parent):
        """计算按钮（第六行）"""
        btn = tk.Button(parent, text="start calculating", command=self._calculate_distance)
        btn.grid(row=6, column=0, pady=10, sticky="ew")

    def _create_txt_controls(self, parent):
        """TXT文件相关控件"""
        # 第七行：新建txt和选择txt
        btn_frame = tk.Frame(parent)
        btn_frame.grid(row=7, column=0, pady=5, sticky="ew")
        tk.Button(btn_frame, text="New txt file", command=self._create_txt_file).pack(side=tk.LEFT, expand=True)
        tk.Button(btn_frame, text="Select txt file", command=self._select_txt_file).pack(side=tk.LEFT, expand=True)

        # 第八行：路径输入框
        path_frame = tk.Frame(parent)
        path_frame.grid(row=8, column=0, pady=5, sticky="ew")

        tk.Label(path_frame, text="txt file path:").pack(side=tk.LEFT)
        entry = tk.Entry(path_frame, textvariable=self.txt_path_var)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # 新增日期输入框
        date_frame = tk.Frame(parent)
        date_frame.grid(row=9, column=0, pady=5, sticky="ew")
        tk.Label(date_frame, text="日期 (YYYYMMDD):").pack(side=tk.LEFT)
        self.date_entry = ttk.Entry(date_frame, width=12)
        self.date_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.date_entry.insert(0, "")  # 初始为空

        # 备注栏
        remark_frame = tk.Frame(parent)
        remark_frame.grid(row=9, column=0, pady=5, sticky="ew")

        tk.Label(remark_frame, text="Remarks:").pack(side=tk.LEFT)
        self.remark_entry = tk.Entry(remark_frame)
        self.remark_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # 写入数据按钮
        write_btn = tk.Button(parent, text="Write to txt file", command=self._write_to_txt)
        write_btn.grid(row=10, column=0, pady=5, sticky="ew")

    def _update_settings_visibility(self):
        """更新设置可见性"""
        current_group = self.group_var.get()
        # 强制更新布局
        self.n_frame.grid(row=4, column=0) if current_group == "group1" else self.n_frame.grid_remove()
        self.s_frame.grid(row=5, column=0) if current_group == "group2" else self.s_frame.grid_remove()

        # 选择器重置逻辑
        if current_group in ["group1", "group2"]:
            self._clear_polygon()
            self._init_selector()
            self.output_label.config(
                text=f"Please start selecting the{'N' if current_group == 'group1' else 'S'}-footpoint area", fg="blue")
        else:
            if hasattr(self, 'selector'):
                self.selector.disconnect_events()

        self.root.update_idletasks()  # 布局刷新

    # ==================================================================
    # ==================文件操作========================
    def select_file(self):  # 选择fits文件并加载数据
        # 重置部分
        if self.data is not None:
            self.data = None
            if self.im:
                self.im.remove()
                self.im = None
            if hasattr(self, 'colorbar') and self.colorbar:
                try:
                    # 安全移除colorbar
                    self.colorbar.remove()
                except (AttributeError, ValueError) as e:
                    print(f"An error occurred while removing the colorbar. {e}")
                # 强制删除引用
                self.colorbar = None
            self.ax.cla()
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.ax.set_frame_on(True)  # 确保axes框架存在

            self.output_label.config(text="")
            self.coordinate_label.config(text="")
            # 重新绘制画布
            self.canvas.draw()
            self._init_selector()

        try:
            filename = askopenfilename(filetypes=[("FITS files", "*.fits")])
            if filename:
                self.current_filename = filename
                self.filename_label.config(text=f"current file：{os.path.basename(filename)}")

                # 使用 with 语句确保文件正确关闭
                with fits.open(filename) as hdul:
                    # 自动查找第一个包含二维数据的HDU
                    for hdu in hdul:
                        if hdu.data is not None and len(hdu.data.shape) == 2:
                            data = hdu.data.astype(np.float32)
                            header = hdu.header
                            break
                    else:
                        raise ValueError("No two-dimensional array data was found in the file")

                    header = self._fix_header(header)

                    self.smap = Map(data, header)

                    # 关键修复：确保太阳半径和观察者坐标存在
                    if 'rsun_ref' not in self.smap.meta:
                        self.smap.meta['rsun_ref'] = 6.957e8
                    if 'rsun_obs' not in self.smap.meta:
                        self.smap.meta['rsun_obs'] = self.smap.meta['rsun_ref']

                    # 确保观察者坐标存在
                    if not hasattr(self.smap, 'observer_coordinate'):
                        from astropy.time import Time
                        self.smap._observer_coordinate = SkyCoord(
                            0 * u.deg, 0 * u.deg,
                            radius=1.495978707e11 * u.m,
                            frame=HeliographicStonyhurst,
                            obstime=Time(self.smap.meta.get('date_obs', self.smap.meta.get('date-obs', '')))
                        )

                    self.data = self.smap.data

                # 填充观测时间
                if hasattr(self.smap, 'date'):
                    obs_datetime = self.smap.date.strftime('%Y-%m-%d %H:%M:%S')
                    self.obs_time_label.config(text=f"observe time：{obs_datetime}")

                    # 自动提取日期并填充到输入框
                    date_str = self.smap.date.strftime('%Y%m%d')
                    self.date_entry.delete(0, tk.END)
                    self.date_entry.insert(0, date_str)

                # 处理 NaN 值
                if np.isnan(self.data).any():
                    min_val = np.nanmin(self.data)
                    self.data = np.nan_to_num(self.data, nan=min_val)

                    # 添加中心点坐标计算
                try:
                    center_y = self.data.shape[0] // 2
                    center_x = self.data.shape[1] // 2

                    if center_y < self.data.shape[0] and center_x < self.data.shape[1]:
                        x_pix = center_x * u.pix
                        y_pix = center_y * u.pix
                        sky_coord = self.smap.pixel_to_world(x_pix, y_pix)

                        if not hasattr(self.smap, 'observer_coordinate'):
                            raise ValueError("FITS header lacks the observer time.")
                        if not hasattr(self.smap, 'date'):
                            raise ValueError("FITS header lacks the observation time.")

                        observer = self.smap.observer_coordinate
                        frame = HeliographicCarrington(obstime=self.smap.date, observer=observer)
                        transformed_coord = sky_coord.transform_to(frame)
                        coord_text = (f"Carrington：\n"
                                      f"longgitude: {transformed_coord.lon:.3f}°\n"
                                      f"latitude: {transformed_coord.lat:.3f}°")
                    else:
                        coord_text = "out of range"
                    self.center_coord_label.config(text=coord_text)
                except Exception as coord_error:
                    self.center_coord_label.config(text=f"calculation error: {str(coord_error)}")

                # 显示图像
                self.update_image()

                """
                # 创建矩形选择器
                self.rect = RectangleSelector(self.ax, self.onselect, useblit=True,
                                          button=[1],  # 只允许左键选择
                                          minspanx=5, minspany=5,
                                          spancoords='pixels',
                                          interactive=True)
                """

                # 多边形选择器
                self.selector = PolygonSelector(
                    self.ax,
                    self.onselect,
                    useblit=True,
                    props=dict(
                        linestyle='-',
                        linewidth=1,
                        color='red',
                        alpha=0.3,
                        marker='o',
                        markersize=5,
                        markerfacecolor='red'
                    )
                )

                self.canvas.mpl_connect('key_press_event', self._handle_key)

                # 连接鼠标点击事件
                self.cid = self.canvas.mpl_connect('button_press_event', self.on_click)


        except Exception as e:
            # 更新 Label 显示错误信息
            self.output_label.config(text=f"发生错误: {str(e)}")

    def _select_txt_file(self):  # 选择txt文件
        file_path = askopenfilename(
            title="Select txt file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if file_path:
            self.txt_path_var.set(file_path)  # 将路径存入输入框
            self.output_label.config(text=f"Selected documents: {file_path}")

    def _create_txt_file(self):  # 创建txt文件
        try:
            desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            file_name = f"distance_results_{timestamp}.txt"
            full_path = os.path.join(desktop_path, file_name)
            with open(full_path, 'w', encoding='utf-8') as f:
                f.write(f"File creation time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                header = "{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<16}\n".format(
                    "DATE", "NPOL", "NARN", "NLON", "NLAT", "SPOL",
                    "SARN", "SLON", "SLAT", "LDEG", "LRAD", "REMARK"
                )
                f.write(header)
            self.txt_path_var.set(full_path)  # 将路径存入输入框
            self.output_label.config(text=f"Created files：\n{full_path}", fg="green")
        except Exception as e:
            self.output_label.config(text=f"Failed to create file：{e}", fg="red")

    def _write_to_txt(self):  # 将计算结果写入txt文件
        target_path = self.txt_path_var.get()  # 从输入框获取路径

        if not target_path:
            self.output_label.config(text="请先创建txt文件", fg="red")
            return

        date_str = self.date_entry.get().strip()
        if not date_str:
            self.output_label.config(text="请输入日期", fg="red")
            return

        try:
            datetime.strptime(date_str, "%Y%m%d")
        except ValueError:
            self.output_label.config(text='please use YYYYMMDD style', fg='red')
            return

        try:
            # 如果文件不存在则创建
            if not os.path.exists(target_path):
                with open(target_path, 'w', encoding='utf-8') as f:
                    f.write(f"creation time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    header = "{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<16}\n".format(
                        "DATE", "NPOL", "NARN", "NLON", "NLAT", "SPOL",
                        "SARN", "SLON", "SLAT", "LDEG", "LRAD", "REMARK"
                    )
                    f.write(header)

            # 获取两组数据
            group1 = self.group_data.get("group1")
            group2 = self.group_data.get("group2")

            # 构建数据行
            line = [date_str]

            # 处理第一组数据
            if group1:
                polarity = "+" if group1.get('polarity') == "positive" else "-"
                lon = f"{group1['carrington'][0]:.2f}"  # 经度
                lat = f"{group1['carrington'][1]:.2f}"  # 纬度
                line.extend([polarity, self.n_ar_entry.get(), lon, lat])
            else:
                line.extend(["-", "-", "-", "-"])

            # 处理第二组数据
            if group2:
                polarity = "+" if group2.get('polarity') == "positive" else "-"
                lon = f"{group2['carrington'][0]:.2f}"
                lat = f"{group2['carrington'][1]:.2f}"
                line.extend([polarity, self.s_ar_entry.get(), lon, lat])
            else:
                line.extend(["-", "-", "-", "-"])

            if hasattr(self, 'distance_deg') and hasattr(self, 'distance_rad'):
                line.extend([f"{self.distance_deg:.3f}", f"{self.distance_rad:.3f}"])
            else:
                line.extend(["-", "-"])

            # 处理备注
            remarks = self.remark_entry.get().strip()

            formatted_line = "{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<8}\t{:<16}\n".format(
                line[0],  # 日期
                line[1],  # N极性
                line[2],  # N_AR
                line[3],  # N经度
                line[4],  # N纬度
                line[5],  # S极性
                line[6],  # S_AR
                line[7],  # S经度
                line[8],  # S纬度
                line[9],  # 距离度
                line[10],  # 距离rad
                remarks  # 备注
            )
            with open(target_path, 'a', encoding='utf-8') as f:
                f.write(formatted_line)

            self.output_label.config(text="wirte success", fg="green")
        except Exception as e:
            self.output_label.config(text=f"write fail: {str(e)}", fg="red")

    # =============================================
    # ==================图像处理=================
    def _fix_header(self, header):
        """修复缺失的太阳图像元数据"""
        # 添加缺失的坐标单位（默认为角秒）
        if 'CUNIT1' not in header:
            header['CUNIT1'] = 'arcsec'
        if 'CUNIT2' not in header:
            header['CUNIT2'] = 'arcsec'

        # 添加坐标类型（Helioprojective Cartesian）
        if 'CTYPE1' not in header:
            header['CTYPE1'] = 'HPLN-TAN'
        if 'CTYPE2' not in header:
            header['CTYPE2'] = 'HPLT-TAN'

        # 添加参考像素位置（图像中心）
        if 'CRPIX1' not in header:
            header['CRPIX1'] = header.get('NAXIS1', 1024) / 2 + 0.5
        if 'CRPIX2' not in header:
            header['CRPIX2'] = header.get('NAXIS2', 1024) / 2 + 0.5

        # 添加参考坐标值（太阳中心）
        if 'CRVAL1' not in header:
            header['CRVAL1'] = 0.0
        if 'CRVAL2' not in header:
            header['CRVAL2'] = 0.0

        # 添加像素尺度（默认值）
        if 'CDELT1' not in header:
            header['CDELT1'] = 1.0  # 弧秒/像素
        if 'CDELT2' not in header:
            header['CDELT2'] = 1.0  # 弧秒/像素

        # 添加太阳半径（标准值）
        if 'RSUN_REF' not in header:
            header['RSUN_REF'] = 6.957e8
        if 'RSUN_OBS' not in header:
            header['RSUN_OBS'] = header['RSUN_REF']

        # 添加观察者位置（地球）
        if 'HGLN_OBS' not in header:
            header['HGLN_OBS'] = 0.0  # 日面经度
        if 'HGLT_OBS' not in header:
            header['HGLT_OBS'] = 0.0  # 日面纬度
        if 'DSUN_OBS' not in header:
            header['DSUN_OBS'] = 1.495978707e11  # 日地距离（米）

        return header

    def update_image(self):  # 更新图像
        try:
            vmin = float(self.vmin_var.get())
            vmax = float(self.vmax_var.get())
            if self.data is not None:
                if self.im is None:
                    # 主图初始化
                    self.ax.cla()
                    self.ax.set_xticks([])
                    self.ax.set_yticks([])

                    # 创建主图图像和颜色条
                    self.im = self.ax.imshow(self.data, vmin=vmin, vmax=vmax, cmap='gray', origin='lower')
                    # self.ax.invert_xaxis() # 翻转x轴
                    self.colorbar = self.fig.colorbar(self.im, ax=self.ax)
                    self.colorbar.ax.set_position([0.85, 0.1, 0.03, 0.8])  # 调整颜色条位置

                    # 设置主图布局
                    self.fig.subplots_adjust(left=0.05, right=0.92)
                else:
                    # 更新主图显示范围
                    self.im.set_clim(vmin, vmax)
                    self.im.norm.vmin = vmin
                    self.im.norm.vmax = vmax

                # 更新颜色条
                if self.colorbar:
                    self.colorbar.update_normal(self.im)

                self.canvas.draw_idle()

        except ValueError:
            self.output_label.config(text="The vmin or vmax input is not a valid number, enter the correct value.")

    def _clear_image(self):  # 清除图像
        """清除所有图像和关联数据"""
        if self.data is not None:
            self.data = None
            if self.im:
                self.im.remove()
                self.im = None
            # 清理colorbar
            if hasattr(self, 'colorbar') and self.colorbar:
                try:
                    self.colorbar.remove()
                except (AttributeError, ValueError):
                    pass
                self.colorbar = None
            # 重置坐标轴（仅保留主图）
            self.ax.cla()

            # 重新创建单图布局
            self.fig.delaxes(self.ax)
            gs = self.fig.add_gridspec(1, 1)  # 改为单图布局
            self.ax = self.fig.add_subplot(gs[0])

            # 调整主图位置
            self.ax.set_position([0.1, 0.1, 0.8, 0.8])
            self.ax.set_xticks([])
            self.ax.set_yticks([])

            # 清除标记
            self._clear_marks()
            # 重绘画布
            self.canvas.draw()
            # 清空数据显示
            self.coordinate_label.config(text="Click coordinates: none")
            self.group1_info.config(text="unselected")
            self.group2_info.config(text="unselected")
            self.distance_label.config(text="uncalculated")

            self.date_entry.delete(0, tk.END)
            self.date_entry.insert(0, 'entry date')
        # 对RectangleSelector的引用进行清理
        if self.rect:
            self.rect.set_visible(False)
            self.rect = None

        # 重置活动区编号
        self.n_ar_entry.delete(0, 'end')
        self.n_ar_entry.insert(0, "0")
        self.s_ar_entry.delete(0, 'end')
        self.s_ar_entry.insert(0, "0")

        # 对事件监听的清理
        if hasattr(self, 'cid'):
            self.canvas.mpl_disconnect(self.cid)

    def _handle_key(self, event):  # 按键处理
        if event.key == 'enter':
            if len(self.selector.verts) > 2:
                self.onselect(self.selector.verts)

            # 清理当前选择器
            self.selector.disconnect_events()
            self.canvas.draw_idle()

            current_group = self.group_var.get()
            if current_group == "group1":
                self.group_var.set("group2")
                self.output_label.config(text="Please continue to select the S footpoint area", fg="blue")
                # 重新初始化选择器
                self._init_selector()
            elif current_group == "group2":
                self.group_var.set("none")
                self.output_label.config(text="The selection is complete and you can start calculating", fg="green")
            self._update_settings_visibility()

    # ============================================
    # ==================坐标计算=================
    def onselect(self, verts):  # 选区
        current_group = self.group_var.get()
        if current_group == "none" or len(verts) < 3:
            return

        try:
            from matplotlib.path import Path
            path = Path(verts)

            # 获取整个数据区域的网格坐标
            ny, nx = self.data.shape
            x, y = np.meshgrid(np.arange(nx), np.arange(ny))
            points = np.vstack((x.ravel(), y.ravel())).T

            # 创建选区掩码
            mask = path.contains_points(points).reshape(ny, nx)
            rows, cols = np.where(mask)

            if len(rows) == 0:
                raise ValueError("The selection does not contain any valid pixels")

            # 清除旧标记
            prefix = 'n' if current_group == 'group1' else 's'
            if hasattr(self, f'marked_point_{current_group}') and getattr(self,
                                                                          f'marked_point_{current_group}') is not None:
                getattr(self, f'marked_point_{current_group}').remove()
                setattr(self, f'marked_point_{current_group}', None)

            # 获取阈值范围
            lower = float(getattr(self, f"{prefix}_lower_threshold_entry").get())
            upper = float(getattr(self, f"{prefix}_upper_threshold_entry").get())

            # 应用阈值过滤
            values = self.data[rows, cols]
            threshold_mask = (values >= lower) & (values <= upper)
            valid_rows = rows[threshold_mask]
            valid_cols = cols[threshold_mask]
            valid_values = values[threshold_mask]

            if len(valid_values) == 0:
                raise ValueError("There is no data in the selection that meets the threshold range")

            # 计算加权平均坐标
            w_x = np.sum(valid_cols * valid_values) / np.sum(valid_values)
            w_y = np.sum(valid_rows * valid_values) / np.sum(valid_values)

            # 转换坐标系
            x_pix = valid_cols * u.pix
            y_pix = valid_rows * u.pix
            sky_coords = self.smap.pixel_to_world(x_pix, y_pix)
            frame_carrington = HeliographicCarrington(obstime=self.smap.date, observer=self.smap.observer_coordinate)
            hc_carrington = sky_coords.transform_to(frame_carrington)

            # 归一化权重
            total_weight = np.sum(valid_values)
            normalized_weights = valid_values / total_weight

            # 计算纬度加权坐标
            w_lat = np.sum(hc_carrington.lat.degree * normalized_weights)

            # 计算经度加权坐标（能够处理跨360度经线情况）
            angles = np.radians(hc_carrington.lon.degree)
            x = np.sum(normalized_weights * np.cos(angles))
            y = np.sum(normalized_weights * np.sin(angles))
            mean_angle = np.arctan2(y, x)
            w_lon = np.degrees(mean_angle) % 360  # 规范化到0-360度

            # 绘制多边形和标记
            self.current_poly = plt.Polygon(verts, closed=True, ec='red', fc='none', lw=1.5)
            self.ax.add_patch(self.current_poly)
            color = 'red' if current_group == 'group1' else 'blue'
            marker = 'x' if current_group == 'group1' else '+'
            """
            setattr(self, f'marked_point_{current_group}', 
                   self.ax.scatter(w_x, w_y, color=color, marker=marker, s=100))
            """
            self.current_marker = self.ax.scatter(w_x, w_y, color=color, marker=marker, s=100)

            # 存储结果
            frame_stony = HeliographicStonyhurst(obstime=self.smap.date)
            hc_stony = sky_coords.transform_to(frame_stony)
            w_lon_stony = np.sum(hc_stony.lon.degree * valid_values) / np.sum(valid_values)
            w_lat_stony = np.sum(hc_stony.lat.degree * valid_values) / np.sum(valid_values)

            self.group_data[current_group] = {
                "pixel": (w_x, w_y),
                "stonyhurst": (w_lon_stony, w_lat_stony),
                "carrington": (w_lon, w_lat),
                "polarity": self.n_polarity_var.get() if current_group == "group1" else self.s_polarity_var.get()
            }

            self._update_group_info()
            self.canvas.draw_idle()

        except Exception as e:
            self.output_label.config(text=f"error: {str(e)}")
            return
        finally:
            self.canvas.draw_idle()

    def on_click(self, event):  # 鼠标点击事件处理
        if event.inaxes == self.ax:
            # 清除上次标记
            self._clear_all_marks()

            # 获取点击坐标
            x, y = event.xdata, event.ydata
            # 获取点击位置的整数坐标
            x_int = int(x)
            y_int = int(y)
            coord_info = ""  # 初始化坐标信息

            # 在点击位置添加标记
            color = 'cyan'
            self.click_marker = self.ax.scatter(
                x, y, s=100,
                facecolors='none',
                edgecolors=color,
                linewidths=1.5
            )
            # 立即重绘
            self.canvas.draw_idle()

            try:
                if self.smap and hasattr(self, 'smap'):
                    # 转换为球面坐标
                    x_pix = x_int * u.pix
                    y_pix = y_int * u.pix
                    sky_coord = self.smap.pixel_to_world(x_pix, y_pix)

                    frame_stony = HeliographicStonyhurst(obstime=self.smap.date)
                    frame_carr = HeliographicCarrington(obstime=self.smap.date,
                                                        observer=self.smap.observer_coordinate)

                    coord_stony = sky_coord.transform_to(frame_stony)
                    coord_carr = sky_coord.transform_to(frame_carr)

                    coord_info = (
                        f"\nStonyhurst: longitude {coord_stony.lon:.3f}° latitude {coord_stony.lat:.3f}°"
                        f"\nCarrington: longitude {coord_carr.lon:.3f}° latitude {coord_carr.lat:.3f}°"
                    )
            except Exception as e:
                coord_info = f"\nCoordinate conversion erro: {str(e)}"
            # 检查坐标是否在数据范围内
            if (self.data is not None
                    and 0 <= y_int < self.data.shape[0]
                    and 0 <= x_int < self.data.shape[1]):
                # 获取该坐标处的数据
                value = self.data[y_int, x_int]
                self.coordinate_label.config(text=f"Click on the location coordinates: (x={x:.1f}, y={y:.1f})\n"
                                                  f"data value: {value:.2f}{coord_info}")
            else:
                self.coordinate_label.config(text=f"Click on the location coordinates: (x={x}, y={y})\n"
                                                  "No valid data")

            from sunpy.coordinates import sun
            carrington_lon = sun.carrington_rotation_number(self.smap.date)

    def _calculate_distance(self):  # 足点距离计算
        # 获取两组数据
        group1 = self.group_data["group1"]
        group2 = self.group_data["group2"]

        if not group1 or not group2:
            self.distance_label.config(text="It is required that both sets of data have been selected!", fg="red")
            return

        try:
            frame = HeliographicCarrington(
                obstime=self.smap.date,
                observer=self.smap.observer_coordinate
            )

            # 创建坐标对象
            coord1 = SkyCoord(
                lon=group1["carrington"][0] * u.deg,
                lat=group1["carrington"][1] * u.deg,
                frame=frame)
            coord2 = SkyCoord(
                lon=group2["carrington"][0] * u.deg,
                lat=group2["carrington"][1] * u.deg,
                frame=frame)

            # 计算角度距离
            sep_deg = coord1.separation(coord2).degree
            # 取最小角度差
            sep_deg = min(sep_deg, 360 - sep_deg)
            # 转换为弧度
            sep_rad = np.deg2rad(sep_deg)

            self.distance_label.config(
                text=f"Spherical angle distance: {sep_deg:.3f}°\narc distance: {sep_rad:.3f} rad",
                fg="black"
            )

            self.distance_deg = sep_deg
            self.distance_rad = sep_rad
        except Exception as e:
            self.distance_label.config(text=f"miscalculation: {str(e)}", fg="red")

    # =================================================
    # ==================数据更新==================
    def _update_group_info(self):  # 更新两组信息显示
        for group in ["group1", "group2"]:
            data = self.group_data[group]
            text = f"{'N-FootPoint' if group == 'group1' else 'S-FootPoint'}: "
            if data:
                ar_number = self.n_ar_entry.get() if group == "group1" else self.s_ar_entry.get()
                text += (f"\nAR-Number:{ar_number}"
                         f"\npixel coordinates: ({data['pixel'][0]:.1f}, {data['pixel'][1]:.1f})"
                         f"\nStonyhurst:\nlongitude: {data['stonyhurst'][0]:.3f}°\nlatitude: {data['stonyhurst'][1]:.3f}°"
                         f"\nCarrington:\nlongitude: {data['carrington'][0]:.3f}°\nlatitude: {data['carrington'][1]:.3f}°")
            else:
                text += "unselected"
            getattr(self, f"{group}_info").config(text=text)

    def _update_polarity(self, lower, upper, prefix):  # 通用极性更新方法
        """通用极性更新方法"""
        getattr(self, f"{prefix}_lower_threshold_entry").delete(0, 'end')
        getattr(self, f"{prefix}_lower_threshold_entry").insert(0, lower)
        getattr(self, f"{prefix}_upper_threshold_entry").delete(0, 'end')
        getattr(self, f"{prefix}_upper_threshold_entry").insert(0, upper)
        self.update_image()

    def _update_polarity_group(self, prefix):  # 通用极性组更新方法
        # 参数校验
        if not hasattr(self, f"{prefix}_polarity_var"):
            raise ValueError(f"invalid prefix：{prefix}")
        polarity_var = getattr(self, f"{prefix}_polarity_var")
        lower, upper = ("10", "10000") if polarity_var.get() == "positive" else ("-10000", "-10")
        self._update_polarity(lower, upper, prefix)

    def _init_selector(self):  # 初始化选择器
        """初始化选择器"""
        if hasattr(self, 'selector'):
            self.selector.disconnect_events()
        self.selector = PolygonSelector(
            self.ax,
            self.onselect,
            useblit=True,
            props=dict(
                linestyle='-',
                linewidth=2,
                color='red',
                alpha=0.5,
                marker='o',
                markersize=5,
                markerfacecolor='red'
            )
        )
        self.canvas.mpl_connect('key_press_event', self._handle_key)

    # ==============================================
    # ====================清理关闭======================
    def _clear_marks(self):  # 清除足点选择标记
        current_group = self.group_var.get()
        if current_group in ["group1", "group2"]:
            attr_name = f"marked_point_{current_group}"
            if hasattr(self, attr_name) and getattr(self, attr_name):
                getattr(self, attr_name).remove()
                setattr(self, attr_name, None)
        self.canvas.draw()

    def _clear_all_marks(self):  # 清除所有标记
        self._clear_marks()
        if hasattr(self, 'click_marker'):
            try:
                self.click_marker.remove()
            except (NotImplementedError, ValueError):
                for i in reversed(range(len(self.ax.collections))):
                    if self.ax.collections[i] == self.click_marker:
                        del self.ax.collections[i]
                        break
            finally:
                if hasattr(self, 'click_marker'):
                    del self.click_marker
        self.canvas.draw_idle()

    def on_close(self):  # 关闭窗口
        # 关闭 matplotlib 图形窗口
        plt.close('all')
        # 关闭窗口时销毁根窗口
        self.root.destroy()

    def _clear_polygon(self):  # 清除多边形选框
        """清除多边形选框"""
        if hasattr(self, 'selector'):
            self.selector.disconnect_events()
            self.selector.clear()
            if hasattr(self, 'poly'):
                self.poly.remove()
                del self.poly
            self.canvas.draw_idle()


# ====================================================
# =====================代码结束=======================

if __name__ == "__main__":
    root = tk.Tk()
    viewer = FITSViewer(root)
    root.mainloop()
