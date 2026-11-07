"""
连接用户已「共享」的 GUI MATLAB 会话，运行 week3_autofocus.m，并把命令窗口输出
记到日记文件，最后回读指标 CSV 做摘要。

前置：用户在 GUI MATLAB 命令窗口执行过一次  matlab.engine.shareEngine
用法：
    python run_autofocus_engine.py              # 自动连最近共享的会话
    python run_autofocus_engine.py MATLAB_96368 # 连指定名称的会话
"""

import os
import sys

import matlab.engine

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DIARY = os.path.join(HERE, "autofocus_run_diary.txt")


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else None
    print("[1/5] 连接共享 MATLAB 会话 ...", flush=True)
    if name:
        eng = matlab.engine.connect_matlab(name)
    else:
        sessions = matlab.engine.find_matlab()
        if not sessions:
            raise SystemExit("未发现共享 MATLAB 会话。请先在 GUI MATLAB 执行: matlab.engine.shareEngine")
        print("    可用会话:", sessions)
        eng = matlab.engine.connect_matlab(sessions[0])
    print("    已连接:", eng, flush=True)

    print("[2/5] 切换工作目录到", HERE, flush=True)
    eng.cd(HERE)

    print("[3/5] 打开日记，运行 week3_autofocus ...", flush=True)
    eng.eval("diary('{}');".format(DIARY.replace("\\", "/")), nargout=0)
    try:
        eng.run("week3_autofocus", nargout=0)
    finally:
        eng.eval("diary off", nargout=0)
    print("    运行结束，日记已写入", DIARY, flush=True)

    print("[4/5] 回读指标 CSV ...", flush=True)
    csv_path = os.path.join(ROOT, "results", "metrics", "autofocus_metrics.csv")
    if os.path.exists(csv_path):
        with open(csv_path, encoding="utf-8") as f:
            print(f.read())
    else:
        print("    (未找到", csv_path, ")")

    print("[5/5] 断开连接（GUI MATLAB 仍在运行）。", flush=True)
    eng.quit()
    print("DONE")


if __name__ == "__main__":
    main()
