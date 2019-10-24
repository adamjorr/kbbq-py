import pytest
import sys
import kbbq.__main__

def test_no_args(monkeypatch, capfd):
    import kbbq.__main__
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]])
        kbbq.main.main()
    captured_noarg = capfd.readouterr()
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0] + "-h"])
        kbbq.main.main()
    captured_h = capfd.readouterr()
    assert captured_noarg == captured_h
