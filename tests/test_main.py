import pytest
import sys
import kbbq.__main__

def test_no_args(monkeypatch, capfd):
    import kbbq.__main__
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]])
        kbbq.__main__.main()
    captured_noarg = capfd.readouterr()
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["-h"])
        with pytest.raises(SystemExit) as e:
            kbbq.__main__.main()
    captured_h = capfd.readouterr()
    assert captured_noarg.out == captured_h.out
