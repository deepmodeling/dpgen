# 配置machine.json将任务提交到lebesgue

您可以登录到lebesgue官方网站 http://lebesgue.dp.tech/ 。后点击左侧的 [功能模块]-[DP-GEN] 后在操作面板中点击左上角的操作手册查看关于DPGEN说明 http://lebesgue.dp.tech/#/jobs/dpgen 。

下面是关于machine.json的配置说明，更多信息或更新请查看官网中的说明。

| field | optional | type | description |
| --- | --- | --- | --- |
| email | false | string | 您的登陆邮箱 |
| password | false | string | 您的lebesgue登陆密码，注意不是服务器的登陆密码 |
| program_id | false | int | 您的项目id(int) 您可以访问 [Lebesgue官网]-[项目管理]-[Program ID]来查看您的项目id |
| job_group_id | true | int | 您在这里配置任务在哪个job下运行。这样任务可以打包在一个job里 |
| rerun | true | int | 非必填,机器回收后是否重跑，0:不重跑(默认), 1：使用原task，2:新建task |
| job_type | false | string | 任务类型，不需要更改 |
| log_file | true | string | 任务日志输出地址，这样任务运行时的输出可以在lebesgue任务预览中看到 |
| command | false | string | 任务运行命令，不需要更改。DPGEN会自动填写 |
| backward_files | false | list[string] | 任务结束后上传的文件，不需要更改。DPGEN会自动填写 |
| job_name | false | string | 任务名，您可以根据任务的类型进行配置 |
| on_demand | true | int | 是否使用ondemand机器默认为0（不开启，只使用Spot机器） |
| image_name | true/false | int | 镜像名称 platform 为ali/aws必填, sugon非必填 |
| disk_size | true/false | int | 数据盘大小 platform 为ali/aws必填，sugon非必填 |
| scass_type | false | string | dp系统机器配置类型 必填,关于scass_type的信息您可以在lebesgue左侧的导航栏内选择[财务分析]-[价格计算器]中通过选择相应的配置即可查看对应的scass_type |
| instance_group_id | true | int | 实例组id，用于多种机型尝试，按顺序尝试,非必填 |
| platform | false | string | aws 亚马逊云服务平台<br />ali 阿里云<br />sugon 曙光云 |
| region | false | string | 计算机器开启地域 |
| zone | false | string | 计算机器开启地区 |
| grouped | false | bool | 是否将同类型的任务在lebesgue网页中打包显示。 |