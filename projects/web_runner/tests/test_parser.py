from pathlib import Path
import pytest
from jurassic_web.parsers import read_table, chart, summary


def test_radiance_table(tmp_path):
    path=tmp_path/'rad.tab'
    path.write_text(''.join(f'# ${i} = Example column {i}\n' for i in range(1,13))+'0 700 0 0 0 0 -8 0 0 0 251 0.5\n')
    table, headers=read_table(path)
    assert table.shape==(1,12)
    assert len(headers)==12
    info=summary(path,['667.7820'])
    assert info['rows']==1
    rad=chart(path,'nadir',['667.7820'])
    assert rad['traces'][0]['y']==[251.0]
    trans=chart(path,'nadir',['667.7820'],'transmittance')
    assert trans['traces'][0]['y']==[0.5]


def test_bad_radiance_rejected(tmp_path):
    path=tmp_path/'rad.tab'
    path.write_text('0 0 nan\n')
    with pytest.raises(ValueError): read_table(path)
