import sys
import os
import imp
import datetime

import numpy as np

from PyQt5.QtWidgets import *
from PyQt5.uic import *

import pyqtgraph as pg

from SeisCore.BinaryFile.BinaryFile import BinaryFile
from SeisCore.Functions.Spectrum import average_spectrum


class Program:
    def __init__(self):
        self.__app = QApplication(sys.argv)
        self.__window = QMainWindow()
        module_folder = imp.find_module('SummationSignal')[1]
        self.__module_folder = module_folder
        ui_path = os.path.join(self.__module_folder,
                               'GUI', 'Forms', 'MainForm.ui')
        self.__ui = loadUi(ui_path, self.__window)
        self.gl = self.__ui.glGraph

        self.channel_a_plot = pg.PlotWidget()
        self.channel_b_plot = pg.PlotWidget()
        self.channel_c_plot = pg.PlotWidget()

        self.channel_a_spectrum = pg.PlotWidget()
        self.channel_b_spectrum = pg.PlotWidget()
        self.channel_c_spectrum = pg.PlotWidget()

        self.gl.addWidget(self.channel_a_plot, 0, 0)
        self.gl.addWidget(self.channel_b_plot, 1, 0)
        self.gl.addWidget(self.channel_c_plot, 2, 0)

        self.gl.addWidget(self.channel_a_spectrum, 0, 1)
        self.gl.addWidget(self.channel_b_spectrum, 1, 1)
        self.gl.addWidget(self.channel_c_spectrum, 2, 1)

        self.__ui.pbFolderExport.clicked.connect(self.select_output_folder)
        self.__ui.pbSelectFiles.clicked.connect(self.select_files)

        self.__ui.pbCalculation.clicked.connect(self.calculation)
        self.__ui.pbSpectrumRefresh.clicked.connect(self.change_spectrum_limits)
        self.__window.show()
        self.__app.exec()

    @property
    def _ui(self):
        return self.__ui

    def select_output_folder(self):
        file_dialog = QFileDialog()
        folder_path = file_dialog.getExistingDirectory()
        self._ui.leOutputFolder.setText(folder_path)

    def select_files(self):
        self._ui.teFilesList.clear()
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        names = file_dialog.getOpenFileNames()[0]
        folder_path = None
        min_dt = None
        max_dt = None
        max_frequency = None
        for path in names:
            if folder_path is None:
                folder_path = os.path.dirname(path)
                self._ui.leInputFolder.setText(folder_path)
            file_name = os.path.basename(path)
            t = file_name.split('.')
            if len(t) != 2:
                continue
            name, extension = t
            if extension not in ('00', 'xx'):
                continue

            bin_data = BinaryFile()
            bin_data.path = path
            bin_data.record_type = 'XYZ'

            is_correct, _ = bin_data.check_correct()
            if not is_correct:
                continue

            if min_dt is None:
                min_dt = bin_data.datetime_start
                max_dt = bin_data.datetime_stop
                max_frequency = bin_data.signal_frequency
            else:
                min_dt = max(min_dt, bin_data.datetime_start)
                max_dt = min(max_dt, bin_data.datetime_stop)
                max_frequency = max(max_frequency, bin_data.signal_frequency)

            self._ui.teFilesList.append(file_name)

        if min_dt is None:
            return None

        if min_dt >= max_dt:
            self._ui.teFilesList.clear()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Отсутствует перекрытие по времени")
            msg.exec_()
        else:
            self._ui.dtMinDt.setMinimumDateTime(min_dt)
            self._ui.dtMinDt.setMaximumDateTime(max_dt)
            self._ui.dtMaxDt.setMinimumDateTime(min_dt)
            self._ui.dtMaxDt.setMaximumDateTime(max_dt)
            self._ui.dtMinDt.setDateTime(min_dt)
            self._ui.dtMaxDt.setDateTime(min_dt + datetime.timedelta(hours=1))
            self._ui.sbResampleFrequency.setMinimum(0)
            self._ui.sbResampleFrequency.setMaximum(max_frequency)
            self._ui.sbResampleFrequency.setValue(max_frequency)

    def calculation(self):
        folder_path = self._ui.leInputFolder.text()
        export_folder_path = self._ui.leOutputFolder.text()

        file_list = self._ui.teFilesList.toPlainText()
        dt_start = self._ui.dtMinDt.dateTime()
        dt_stop = self._ui.dtMaxDt.dateTime()
        if folder_path is None or file_list is None or export_folder_path is None:
            return None
        files = file_list.split('\n')

        extract_data = None
        for index, file in enumerate(files):
            full_path = os.path.join(folder_path, file)

            bin_data = BinaryFile()
            bin_data.path = full_path
            bin_data.record_type = 'XYZ'
            bin_data.use_avg_values = False
            bin_data.read_date_time_start = dt_start.toPyDateTime()
            bin_data.read_date_time_stop = dt_stop.toPyDateTime()
            bin_data.resample_frequency = self._ui.sbResampleFrequency.value()
            if bin_data.signals is None:
                continue

            if extract_data is None:
                discrete_count = bin_data.signals.shape[0]
                extract_data = np.empty(shape=(discrete_count, 3 * len(files)),
                                        dtype=np.int)
            extract_data[:, 3 * index:3 * (index + 1)] = bin_data.signals

        if extract_data is None:
            return None

        sum_trace = np.empty(shape=(extract_data.shape[0], 3),
                             dtype=np.int)
        channel_plots = [self.channel_a_plot, self.channel_b_plot, self.channel_c_plot]
        spectrum_plots = [self.channel_a_spectrum, self.channel_b_spectrum,
                          self.channel_c_spectrum]

        window_size = self._ui.sbWindowSize.value()
        overlap_size = self._ui.sbOverlapSize.value()
        median_filter = self._ui.sbMedianFilter.value()
        marmett_filter = self._ui.sbMarmettFilter.value()
        resample_frequency = self._ui.sbResampleFrequency.value()
        if marmett_filter == 0:
            marmett_filter = None
        if median_filter == 0:
            median_filter = None
        spectrums = None

        for i in range(3):
            channel_signals = extract_data[:, i:3 * len(files):3]
            sum_trace[:, i] = np.sum(channel_signals, axis=1)
            channel_plot = channel_plots[i]
            spectrum_plot = spectrum_plots[i]
            channel_plot.clear()
            channel_plot.plot(sum_trace[:, i], pen=(196, 255, 18))
            spectrum_plot.clear()

            spectrum_data = average_spectrum(signal=sum_trace[:, i],
                                             frequency=resample_frequency,
                                             window=window_size,
                                             overlap=overlap_size,
                                             med_filter=median_filter,
                                             marmett_filter=marmett_filter)
            if spectrums is None:
                spectrums = np.empty(shape=(spectrum_data.shape[0], 4),
                                     dtype=np.float)
                spectrums[:, :2] = spectrum_data
            else:
                spectrums[:, 1 + i] = spectrum_data[:, 1]

            spectrum_plot.plot(spectrum_data, pen=(255, 0, 0))

        tmp_file_a = os.path.join(export_folder_path, 'sum_trace.npy')
        tmp_file_b = os.path.join(export_folder_path, 'spectrums.npy')
        np.save(tmp_file_a, sum_trace)
        np.save(tmp_file_b, spectrums)

        output_file_sum_trace = os.path.join(export_folder_path, 'SumTraces.dat')
        output_file_spectrums = os.path.join(export_folder_path, 'Spectrums.dat')
        np.savetxt(output_file_sum_trace, sum_trace, fmt='%i', delimiter='\t',
                   header='Channel_1\tChannel_2\tChannel_3', comments='')
        np.savetxt(output_file_spectrums, spectrums, fmt='%f', delimiter='\t',
                   header='Frequency\tChannel_1\tChannel_2\tChannel_3',
                   comments='')

    def change_spectrum_limits(self):
        export_folder_path = self._ui.leOutputFolder.text()
        tmp_file_b = os.path.join(export_folder_path, 'spectrums.npy')

        if not os.path.exists(tmp_file_b):
            return None
        min_freq = self._ui.deMinFrequency.value()
        max_freq = self._ui.deMaxFrequency.value()

        spec_data = np.load(tmp_file_b)
        spec_data = spec_data[(spec_data[:, 0] >= min_freq) &
                              (spec_data[:, 0] <= max_freq)]
        spectrum_plots = [self.channel_a_spectrum, self.channel_b_spectrum,
                          self.channel_c_spectrum]
        for i in range(3):
            spectrum_plot = spectrum_plots[i]
            spectrum_plot.clear()
            selection = np.empty(shape=(spec_data.shape[0], 2))
            selection[:, 0] = spec_data[:, 0]
            selection[:, 1] = spec_data[:, 1 + i]
            spectrum_plot.plot(selection, pen=(255, 0, 0))
            spectrum_plot.setRange(
                xRange=(selection[0, 0], selection[-1, 0]),
                yRange=(min(selection[:, 1]), max(selection[:, 1])))


def run():
    Program()

