import tkinter as tk
from tkinter import ttk, messagebox
import serial
import serial.tools.list_ports
import time

BAUDRATE = 9600
IC_COUNT = 8
BITS_PER_IC = 4


class UartGui:
    def __init__(self, root):
        self.root = root
        self.root.title("UART IC Controller")
        self.root.geometry("980x700")
        self.root.minsize(980, 700)

        self.ser = None
        self.port_map = {}
        self.rx_buffer = ""

        self.port_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Bağlı değil")
        self.frame_var = tk.StringVar()

        self.ic_vars = [tk.StringVar(value="0000") for _ in range(IC_COUNT)]

        self.build_ui()
        self.refresh_ports()
        self.update_frame_preview()

        # Arka planda sürekli UART oku
        self.root.after(100, self.poll_serial)

    def build_ui(self):
        main = ttk.Frame(self.root, padding=10)
        main.pack(fill="both", expand=True)

        conn_frame = ttk.LabelFrame(main, text="Bağlantı", padding=10)
        conn_frame.pack(fill="x", pady=(0, 10))

        ttk.Label(conn_frame, text="COM Port:").grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.port_combo = ttk.Combobox(
            conn_frame,
            textvariable=self.port_var,
            state="readonly",
            width=45
        )
        self.port_combo.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        ttk.Button(conn_frame, text="Portları Yenile", command=self.refresh_ports).grid(
            row=0, column=2, padx=5, pady=5
        )
        ttk.Button(conn_frame, text="Bağlan", command=self.connect_port).grid(
            row=0, column=3, padx=5, pady=5
        )
        ttk.Button(conn_frame, text="Bağlantıyı Kes", command=self.disconnect_port).grid(
            row=0, column=4, padx=5, pady=5
        )

        ttk.Label(conn_frame, text="Durum:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        ttk.Label(conn_frame, textvariable=self.status_var).grid(
            row=1, column=1, columnspan=4, padx=5, pady=5, sticky="w"
        )

        ic_frame = ttk.LabelFrame(main, text="Entegre Değerleri", padding=10)
        ic_frame.pack(fill="x", pady=(0, 10))

        ttk.Label(ic_frame, text="Entegre", width=12).grid(row=0, column=0, padx=5, pady=5)
        ttk.Label(ic_frame, text="4 Bit", width=12).grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(ic_frame, text="Hızlı", width=20).grid(row=0, column=2, padx=5, pady=5)

        vcmd = (self.root.register(self.validate_bits), "%P")

        for i in range(IC_COUNT):
            ttk.Label(ic_frame, text=f"IC {i+1}").grid(row=i + 1, column=0, padx=5, pady=6, sticky="w")

            entry = ttk.Entry(
                ic_frame,
                textvariable=self.ic_vars[i],
                width=10,
                justify="center",
                validate="key",
                validatecommand=vcmd
            )
            entry.grid(row=i + 1, column=1, padx=5, pady=6, sticky="w")
            entry.bind("<KeyRelease>", lambda e: self.update_frame_preview())

            quick_frame = ttk.Frame(ic_frame)
            quick_frame.grid(row=i + 1, column=2, padx=5, pady=6, sticky="w")

            ttk.Button(
                quick_frame,
                text="0000",
                command=lambda idx=i: self.set_ic(idx, "0000")
            ).pack(side="left", padx=2)

            ttk.Button(
                quick_frame,
                text="1111",
                command=lambda idx=i: self.set_ic(idx, "1111")
            ).pack(side="left", padx=2)

        preview_frame = ttk.LabelFrame(main, text="Oluşan UART Frame", padding=10)
        preview_frame.pack(fill="x", pady=(0, 10))

        preview_entry = ttk.Entry(preview_frame, textvariable=self.frame_var, font=("Consolas", 12))
        preview_entry.pack(fill="x", padx=5, pady=5)

        button_frame = ttk.Frame(main)
        button_frame.pack(fill="x", pady=(0, 10))

        ttk.Button(button_frame, text="Hepsini 0 Yap", command=self.set_all_zero).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Hepsini 1 Yap", command=self.set_all_one).pack(side="left", padx=5)
        ttk.Button(button_frame, text="PING Gönder", command=self.send_ping).pack(side="left", padx=5)
        ttk.Button(button_frame, text="HELLO Gönder", command=self.send_hello).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Test Frame Gönder", command=self.send_test_frame).pack(side="left", padx=5)

        self.send_button = tk.Button(
            button_frame,
            text="GÖNDER",
            font=("Arial", 14, "bold"),
            height=2,
            width=16,
            command=self.send_frame
        )
        self.send_button.pack(side="right", padx=5)

        log_frame = ttk.LabelFrame(main, text="Log", padding=10)
        log_frame.pack(fill="both", expand=True)

        self.log_text = tk.Text(log_frame, height=16, font=("Consolas", 10))
        self.log_text.pack(fill="both", expand=True)

    def log(self, msg):
        stamp = time.strftime("%H:%M:%S")
        self.log_text.insert("end", f"[{stamp}] {msg}\n")
        self.log_text.see("end")

    def validate_bits(self, new_value):
        if new_value == "":
            return True
        if len(new_value) > BITS_PER_IC:
            return False
        return all(ch in "01" for ch in new_value)

    def normalize_bits(self, value):
        value = value.strip()
        if value == "":
            return "0000"
        if not all(ch in "01" for ch in value):
            raise ValueError("Sadece 0 ve 1 girilebilir.")
        if len(value) > 4:
            raise ValueError("Her IC için en fazla 4 bit girilebilir.")
        return value.zfill(4)

    def set_ic(self, idx, value):
        self.ic_vars[idx].set(value)
        self.update_frame_preview()

    def set_all_zero(self):
        for var in self.ic_vars:
            var.set("0000")
        self.update_frame_preview()

    def set_all_one(self):
        for var in self.ic_vars:
            var.set("1111")
        self.update_frame_preview()

    def build_frame(self):
        values = []
        for i, var in enumerate(self.ic_vars):
            try:
                bits = self.normalize_bits(var.get())
            except ValueError as e:
                raise ValueError(f"IC {i+1}: {e}")
            values.append(bits)
        return "[" + ",".join(values) + "]"

    def update_frame_preview(self):
        try:
            self.frame_var.set(self.build_frame())
        except ValueError:
            pass

    def refresh_ports(self):
        self.port_map.clear()
        port_labels = []

        for p in serial.tools.list_ports.comports():
            desc = p.description if p.description else "Açıklama yok"
            label = f"{p.device} - {desc}"
            port_labels.append(label)
            self.port_map[label] = p.device

        self.port_combo["values"] = port_labels

        if port_labels:
            if self.port_var.get() not in port_labels:
                self.port_var.set(port_labels[0])
        else:
            self.port_var.set("")

        self.log(f"Bulunan portlar: {port_labels if port_labels else 'Yok'}")

    def get_selected_device(self):
        selected = self.port_var.get().strip()
        if not selected:
            return ""
        return self.port_map.get(selected, selected.split(" - ")[0].strip())

    def connect_port(self):
        port = self.get_selected_device()
        if not port:
            messagebox.showerror("Hata", "Önce bir COM port seç.")
            return

        try:
            if self.ser and self.ser.is_open:
                self.ser.close()

            self.ser = serial.Serial(
                port=port,
                baudrate=BAUDRATE,
                bytesize=serial.EIGHTBITS,
                parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE,
                timeout=0.05,
                write_timeout=1
            )

            time.sleep(0.3)

            self.status_var.set(f"Bağlı: {port} @ {BAUDRATE}")
            self.log(f"Bağlandı -> {port} @ {BAUDRATE}")

        except Exception as e:
            self.status_var.set("Bağlı değil")
            messagebox.showerror("Bağlantı hatası", str(e))
            self.log(f"Bağlantı hatası: {e}")

    def disconnect_port(self):
        try:
            if self.ser and self.ser.is_open:
                port_name = self.ser.port
                self.ser.close()
                self.log(f"Bağlantı kesildi -> {port_name}")
        except Exception as e:
            self.log(f"Bağlantı kesme hatası: {e}")

        self.status_var.set("Bağlı değil")

    def poll_serial(self):
        try:
            if self.ser and self.ser.is_open:
                waiting = self.ser.in_waiting
                if waiting > 0:
                    data = self.ser.read(waiting).decode("ascii", errors="ignore")
                    self.rx_buffer += data

                    while "\n" in self.rx_buffer:
                        line, self.rx_buffer = self.rx_buffer.split("\n", 1)
                        line = line.strip("\r").strip()
                        if line:
                            self.log(f"RX: {line}")
        except Exception as e:
            self.log(f"UART okuma hatası: {e}")

        self.root.after(100, self.poll_serial)

    def send_raw_text(self, text):
        if not self.ser or not self.ser.is_open:
            messagebox.showerror("Hata", "Önce porta bağlan.")
            return

        try:
            packet = text + "\r\n"
            self.ser.write(packet.encode("ascii"))
            self.ser.flush()
            self.log(f"TX: {text}")
        except Exception as e:
            messagebox.showerror("Gönderim hatası", str(e))
            self.log(f"Gönderim hatası: {e}")

    def send_frame(self):
        try:
            frame = self.build_frame()
        except ValueError as e:
            messagebox.showerror("Geçersiz veri", str(e))
            return

        # MCU parser'ı SET:[....] bekliyor
        self.send_raw_text("SET:" + frame)

    def send_test_frame(self):
        test_frame = "[0000,0000,0000,0000,0000,0000,0000,0000]"
        self.send_raw_text("SET:" + test_frame)

    def send_ping(self):
        self.send_raw_text("PING")

    def send_hello(self):
        self.send_raw_text("HELLO")


def main():
    root = tk.Tk()
    app = UartGui(root)
    root.mainloop()


if __name__ == "__main__":
    main()