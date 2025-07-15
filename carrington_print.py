import tkinter as tk
from tkinter import ttk, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
import numpy as np
from astropy.io import fits
import sunpy.map
import astropy.units as u
from skimage import measure  # 用于绘制选区轮廓
from sunpy.coordinates import frames
import sunpy
class MagMapAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Solar Magnetogram Analyzer (Dual Selection)")

        # 初始化变量
        self.data = None
        self.smap = None
        self.current_vmin = -1000
        self.current_vmax = 1000
        self.positive_pixels = None  # 正极选区
        self.negative_pixels = None  # 负极选区
        self.positive_threshold = 50   # 正极阈值
        self.negative_threshold = -50   # 负极阈值
        self.pos_coord = None  # 正极日面坐标
        self.neg_coord = None  # 负极日面坐标
        # 新增实例变量
        self.observation_date = ""
        self.pos_pixel = (np.nan, np.nan)
        self.neg_pixel = (np.nan, np.nan)
        self.pos_carrington = None
        self.neg_carrington = None
        self.distance_rad = 0.0

        # 创建界面
        self.create_widgets()

        # 初始化Lasso选择器
        self.lasso = LassoSelector(self.ax, self.on_select, useblit=True)
        self.selection_mode = "positive"  # 当前选区模式


    def create_widgets(self):
        # 文件选择部分
        file_frame = ttk.Frame(self.root)
        file_frame.pack(fill=tk.X, padx=5, pady=5)

        self.btn_browse = ttk.Button(file_frame, text="Browse", command=self.load_file)
        self.btn_browse.pack(side=tk.LEFT, padx=5)

        self.lbl_file = ttk.Label(file_frame, text="No file selected")
        self.lbl_file.pack(side=tk.LEFT, padx=5)

        # 模式选择按钮
        mode_frame = ttk.Frame(self.root)
        mode_frame.pack(fill=tk.X, padx=5, pady=2)

        self.btn_positive = ttk.Button(
            mode_frame, text="Select Positive Area",
            command=lambda: self.set_selection_mode("positive")
        )
        self.btn_positive.pack(side=tk.LEFT, padx=5)

        self.btn_negative = ttk.Button(
            mode_frame, text="Select Negative Area",
            command=lambda: self.set_selection_mode("negative")
        )
        self.btn_negative.pack(side=tk.LEFT, padx=5)

        self.btn_clear = ttk.Button(
            mode_frame, text="Clear All",
            command=self.clear_selections
        )
        self.btn_clear.pack(side=tk.LEFT, padx=5)
        # 在 create_widgets 方法的结果展示部分添加：

        # 参数设置
        param_frame = ttk.Frame(self.root)
        param_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(param_frame, text="Vmin:").grid(row=0, column=0, padx=2)
        self.entry_vmin = ttk.Entry(param_frame, width=8)
        self.entry_vmin.insert(0, str(self.current_vmin))
        self.entry_vmin.grid(row=0, column=1, padx=2)

        ttk.Label(param_frame, text="Vmax:").grid(row=0, column=2, padx=2)
        self.entry_vmax = ttk.Entry(param_frame, width=8)
        self.entry_vmax.insert(0, str(self.current_vmax))
        self.entry_vmax.grid(row=0, column=3, padx=2)
        # 添加阈值输入
        ttk.Label(param_frame, text="+Th:").grid(row=0, column=4, padx=2)
        self.entry_pos_th = ttk.Entry(param_frame, width=8)
        self.entry_pos_th.insert(0, str(self.positive_threshold))
        self.entry_pos_th.grid(row=0, column=5, padx=2)

        ttk.Label(param_frame, text="-Th:").grid(row=0, column=6, padx=2)
        self.entry_neg_th = ttk.Entry(param_frame, width=8)
        self.entry_neg_th.insert(0, str(self.negative_threshold))
        self.entry_neg_th.grid(row=0, column=7, padx=2)

        self.btn_apply = ttk.Button(param_frame, text="Apply", command=self.apply_parameters)
        self.btn_apply.grid(row=0, column=8, padx=5)
        self.btn_apply = ttk.Button(param_frame, text="Apply", command=self.apply_parameters)
        self.btn_apply.grid(row=0, column=4, padx=5)

        # 图像显示
        self.fig = plt.figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # 结果展示
        result_frame = ttk.Frame(self.root)
        result_frame.pack(fill=tk.BOTH, padx=5, pady=5)

        # 调整行号以添加新坐标显示
        ttk.Label(result_frame, text="Positive (Pixel):").grid(row=0, column=0, padx=2)
        self.lbl_pos_pixel = ttk.Label(result_frame, text="")
        self.lbl_pos_pixel.grid(row=0, column=1, padx=2)

        ttk.Label(result_frame, text="Negative (Pixel):").grid(row=0, column=2, padx=2)
        self.lbl_neg_pixel = ttk.Label(result_frame, text="")
        self.lbl_neg_pixel.grid(row=0, column=3, padx=2)

        ttk.Label(result_frame, text="Positive (Helio):").grid(row=1, column=0, padx=2)
        self.lbl_pos_helio = ttk.Label(result_frame, text="")
        self.lbl_pos_helio.grid(row=1, column=1, padx=2)

        ttk.Label(result_frame, text="Negative (Helio):").grid(row=1, column=2, padx=2)
        self.lbl_neg_helio = ttk.Label(result_frame, text="")
        self.lbl_neg_helio.grid(row=1, column=3, padx=2)

        # 新增卡林顿坐标显示
        ttk.Label(result_frame, text="Positive (Carrington):").grid(row=2, column=0, padx=2)
        self.lbl_pos_carrington = ttk.Label(result_frame, text="")
        self.lbl_pos_carrington.grid(row=2, column=1, padx=2)

        ttk.Label(result_frame, text="Negative (Carrington):").grid(row=2, column=2, padx=2)
        self.lbl_neg_carrington = ttk.Label(result_frame, text="")
        self.lbl_neg_carrington.grid(row=2, column=3, padx=2)

        ttk.Label(result_frame, text="Distance:").grid(row=3, column=0, padx=2)
        self.lbl_distance = ttk.Label(result_frame, text="")
        self.lbl_distance.grid(row=3, column=1, columnspan=3, padx=2, sticky=tk.W)
        # 在 create_widgets 方法的结果展示部分添加：
        ttk.Button(result_frame, text="Print", command=self.print_json).grid(row=4, column=0, columnspan=4, pady=5)

    def set_selection_mode(self, mode):
        """设置选区模式并更新按钮状态"""
        self.selection_mode = mode
        self.btn_positive.state(["!pressed" if mode == "positive" else "pressed"])
        self.btn_negative.state(["!pressed" if mode == "negative" else "pressed"])

    def clear_selections(self):
        """清除所有选区"""
        self.positive_pixels = None
        self.negative_pixels = None
        self.update_display(clear_contours=True)
        self.clear_results()

    def clear_results(self):
        """清除结果显示"""
        self.lbl_pos_pixel.config(text="")
        self.lbl_neg_pixel.config(text="")
        self.lbl_pos_helio.config(text="")
        self.lbl_neg_helio.config(text="")
        self.lbl_distance.config(text="")

    def load_file(self):
        """加载FITS文件"""
        filepath = filedialog.askopenfilename(filetypes=[("FITS files", "*.fits")])
        if not filepath:
            return

        try:
            with fits.open(filepath) as hdul:
                hdul.verify('fix')  # 处理所有无效关键字（如 DATAKURT）
                # 处理坐标系标签
                header = hdul[0].header.copy()

                if 'CTYPE1' in header:
                    header['CTYPE1'] = 'HPLN-TAN' if 'solar-x' in header['CTYPE1'].lower() else header['CTYPE1']
                if 'CTYPE2' in header:
                     header['CTYPE2'] = 'HPLT-TAN' if 'solar-y' in header['CTYPE2'].lower() else header['CTYPE2']
                # header['DSUN_OBS'] = 1.496e11  # SOHO卫星到太阳的距离（单位：米）
                # header['HGLN_OBS'] = 0.0  # 日面经度（与地球对齐）
                # header['HGLT_OBS'] = 0.0  # 日面纬度（位于黄道面）
                self.smap = sunpy.map.Map(hdul[0].data, header)
                self.data = self.smap.data

            self.lbl_file.config(text=filepath.split("/")[-1])
            self.current_vmin = np.nanpercentile(self.data, 1)
            self.current_vmax = np.nanpercentile(self.data, 99)
            self.entry_vmin.delete(0, tk.END)
            self.entry_vmax.delete(0, tk.END)
            self.entry_vmin.insert(0, f"{self.current_vmin:.2f}")
            self.entry_vmax.insert(0, f"{self.current_vmax:.2f}")
            self.update_display(clear_contours=True)
            self.clear_selections()
            # 替换原有日期打印代码
            self.observation_date = header.get('DATE-OBS', header.get('DATE-BEG', ''))
            print("FITS Header:", self.smap.meta)
            date_obs = header.get('DATE-OBS', header.get('DATE-BEG', 'UNKNOWN'))
            print(f"观测日期: {date_obs}")
        except Exception as e:
            print(f"Error loading file: {str(e)}")

    def apply_parameters(self):
        """应用显示范围参数"""
        try:
            self.current_vmin = float(self.entry_vmin.get())
            self.current_vmax = float(self.entry_vmax.get())
            self.positive_threshold = float(self.entry_pos_th.get())
            self.negative_threshold = float(self.entry_neg_th.get())
            self.update_display()
        except ValueError:
            pass

    def update_display(self, clear_contours=False):
        """更新图像显示"""
        self.ax.clear()

        if self.data is not None:
            self.ax.imshow(self.data, cmap='gray', origin='lower',
                           vmin=self.current_vmin, vmax=self.current_vmax)

            # 绘制选区轮廓
            if not clear_contours:
                self.draw_contours(self.positive_pixels, 'red')
                self.draw_contours(self.negative_pixels, 'blue')

        self.canvas.draw()

    def draw_contours(self, mask, color):
        """绘制选区轮廓"""
        if mask is None:
            return
        try:
            contours = measure.find_contours(mask.astype(int), 0.5)
            for contour in contours:
                self.ax.plot(contour[:, 1], contour[:, 0], linewidth=1.5,
                             color=color, linestyle='--')
        except Exception as e:
            print(f"Error drawing contour: {str(e)}")

    def on_select(self, verts):
        """处理Lasso选区"""
        if self.data is None:
            return

        path = Path(verts)
        h, w = self.data.shape
        x, y = np.meshgrid(np.arange(w), np.arange(h))
        mask = path.contains_points(np.vstack((x.ravel(), y.ravel())).T).reshape(h, w)

        # 根据模式保存选区
        if self.selection_mode == "positive":
            self.positive_pixels = mask
        else:
            self.negative_pixels = mask

        self.update_display()
        self.calculate_centroids()

    def calculate_centroids(self):
        """更新后的质心计算入口"""
        # 计算正极质心
        pos_x, pos_y = self.calculate_single_centroid(self.positive_pixels, True)
        self.lbl_pos_pixel.config(text=f"({pos_x:.1f}, {pos_y:.1f})" if not np.isnan(pos_x) else "N/A")

        # 计算负极质心
        neg_x, neg_y = self.calculate_single_centroid(self.negative_pixels, False)
        self.lbl_neg_pixel.config(text=f"({neg_x:.1f}, {neg_y:.1f})" if not np.isnan(neg_x) else "N/A")

        # 转换坐标系
        pos_helio, self.pos_coord, pos_carrington, _ = self.to_helio(pos_x, pos_y)
        self.lbl_pos_helio.config(text=pos_helio)
        self.lbl_pos_carrington.config(text=pos_carrington)

        neg_helio, self.neg_coord, neg_carrington, _ = self.to_helio(neg_x, neg_y)
        self.lbl_neg_helio.config(text=neg_helio)
        self.lbl_neg_carrington.config(text=neg_carrington)

        # 计算角度距离
        if self.pos_coord and self.neg_coord:
            theta_deg, theta_rad = self.calculate_angular_distance(
                self.pos_coord, self.neg_coord
            )
            self.lbl_distance.config(
                text=f"{theta_deg:.2f}° / {theta_rad:.4f} rad"
            )
        else:
            self.lbl_distance.config(text="N/A")

        # 替换原有坐标显示代码
        # 正极部分
        self.pos_pixel = (pos_x, pos_y)
        pos_helio, self.pos_coord, pos_carrington, self.pos_carrington = self.to_helio(pos_x, pos_y)

        # 负极部分
        self.neg_pixel = (neg_x, neg_y)
        neg_helio, self.neg_coord, neg_carrington, self.neg_carrington = self.to_helio(neg_x, neg_y)

        # 计算距离
        if self.pos_coord and self.neg_coord:
            self.distance_rad = self.calculate_angular_distance(self.pos_coord, self.neg_coord)[1]
        else:
            self.distance_rad = 0.0

    def calculate_single_centroid(self, mask, is_positive):
        """根据选区和阈值计算质心"""
        if mask is None:
            return np.nan, np.nan

        # 获取选区内的像素坐标
        rows, cols = np.where(mask)
        if len(rows) == 0:
            return np.nan, np.nan

        # 提取磁通量数据并应用阈值过滤
        data_values = self.data[rows, cols]
        if is_positive:
            valid_mask = data_values >= self.positive_threshold
        else:
            valid_mask = data_values <= self.negative_threshold

        valid_rows = rows[valid_mask]
        valid_cols = cols[valid_mask]

        if len(valid_rows) == 0:
            return np.nan, np.nan

        # 计算加权质心
        weights = np.abs(self.data[valid_rows, valid_cols])
        total_weight = np.sum(weights)

        if total_weight == 0:
            return np.nan, np.nan

        return (np.sum(valid_cols * weights) / total_weight,
                np.sum(valid_rows * weights) / total_weight)

    def to_helio(self, x_pixel, y_pixel):
        """转换像素坐标到日面坐标和卡林顿坐标"""
        if np.isnan(x_pixel) or np.isnan(y_pixel):
            return "N/A", None, "N/A", None

        try:
            coord = self.smap.pixel_to_world(x_pixel * u.pixel, y_pixel * u.pixel)

            # 转换为斯通赫斯特坐标系
            stonyhurst_coord = coord.transform_to(frames.HeliographicStonyhurst)
            lon_stonyhurst = stonyhurst_coord.lon.wrap_at(180 * u.deg).to(u.deg).value
            lat_stonyhurst = stonyhurst_coord.lat.to(u.deg).value

            # 转换为卡林顿坐标系
            carrington_coord = coord.transform_to(frames.HeliographicCarrington)
            lon_carrington = carrington_coord.lon.wrap_at(360 * u.deg).to(u.deg).value
            lat_carrington = carrington_coord.lat.to(u.deg).value

            # 检查是否在日面范围内
            if np.sqrt(coord.Tx  **  2 + coord.Ty  **  2) > self.smap.rsun_obs:
                return "Out of disk", None, "Out of disk", None

            stonyhurst_str = f"{lon_stonyhurst:.2f}°, {lat_stonyhurst:.2f}°"
            carrington_str = f"{lon_carrington:.2f}°, {lat_carrington:.2f}°"
            return stonyhurst_str, (lon_stonyhurst, lat_stonyhurst), carrington_str, (lon_carrington, lat_carrington)
        except Exception as e:
            print(f"Coordinate conversion error: {str(e)}")
            return "N/A", None, "N/A", None

    def calculate_angular_distance(self, coord1, coord2):
        """计算两个日面坐标之间的中心角"""
        if coord1 is None or coord2 is None:
            return (np.nan, np.nan)

        try:
            # 将经纬度转换为弧度
            lon1, lat1 = np.deg2rad(coord1)
            lon2, lat2 = np.deg2rad(coord2)

            # 转换为三维笛卡尔坐标（单位球面）
            x1 = np.cos(lat1) * np.cos(lon1)
            y1 = np.cos(lat1) * np.sin(lon1)
            z1 = np.sin(lat1)

            x2 = np.cos(lat2) * np.cos(lon2)
            y2 = np.cos(lat2) * np.sin(lon2)
            z2 = np.sin(lat2)

            # 计算点积并确定中心角
            dot_product = x1 * x2 + y1 * y2 + z1 * z2
            theta_rad = np.arccos(np.clip(dot_product, -1, 1))
            theta_deg = np.rad2deg(theta_rad)

            return theta_deg, theta_rad
        except Exception as e:
            print(f"Distance calculation error: {str(e)}")
            return (np.nan, np.nan)

    def print_json(self):
        import json

        def format_coord(coord):
            return [round(coord[0], 2), round(coord[1], 2)] if coord else []

        output = {
            "date": self.observation_date.split("T")[0] if self.observation_date else "",
            "positive": {
                "pixel_centroid": [round(self.pos_pixel[0], 1), round(self.pos_pixel[1], 1)]
                if not np.isnan(self.pos_pixel[0]) else [],
                "stonyhurst": format_coord(self.pos_coord),
                "carrington": format_coord(self.pos_carrington)
            },
            "negative": {
                "pixel_centroid": [round(self.neg_pixel[0], 1), round(self.neg_pixel[1], 1)]
                if not np.isnan(self.neg_pixel[0]) else [],
                "stonyhurst": format_coord(self.neg_coord),
                "carrington": format_coord(self.neg_carrington)
            },
            "distance_rad": round(self.distance_rad, 4)
        }

        print(json.dumps(output, indent=4))
if __name__ == "__main__":

    root = tk.Tk()
    app = MagMapAnalyzer(root)
    root.mainloop()
