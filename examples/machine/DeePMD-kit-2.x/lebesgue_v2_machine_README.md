# Config machine.json file in order to submit task to lebesgue platform.

You can login to lebesgue official website http://lebesgue.dp.tech/ . Then click [Function]-[DPGEN]-[manual]\(On the top left corner of the function panel\) from left navigator bar http://lebesgue.dp.tech/#/jobs/dpgen. 

Below is the description of each json fields, please visit official documentation for more information and update.

| field | optional | type | description |
| --- | --- | --- | --- |
| email | false | string | your lebesgue login email |
| password | false | string | your lebesgue login password (note this is not your remote machine login password) |
| program_id | false | int | your program id(int) you can find it in [Lebesgue Official website]-[Programs]-[Program ID] to view your program id |
| api_version | true | int| (default 1) the api_version inside input_data is different from the outside one, which is used to decide which api version will be called to lebesgue. lebesgue currently support version 1 and 2, and version 1 will be deprecate in the future. |
| job_group_id | true | int | config this to specific job_group so submitted jobs can be view as a whole group in the webpage.
| rerun | true | int | if the submitted job terminate unsuccessfully, does it need to be rerun.
| job_type | false | string | job type, should be indicate |
| log_file | true | string | the location of log file, where you can view the log in webpage |
| job_name | false | string | job group name  |
| on_demand | true | int | default:0, 0:use spot machine 1:use ondemand machine |
| image_name | true/false | int | image name, necessary when platform is ali or aws, optional when platform is sugon |
| disk_size | true/false | int | disk size (GB), necessary when platform is ali or aws, optional when platform is sugon |
| scass_type | false | string | machine configuration, about scass_type, you can find them on [lebesgue official website] - [Finance]-[Price calculator] to select disire machine configuration. invalid when instance_group_id is present |
| instance_group_id | true | int | group of scass type |
| platform | false | string | avaliable platform: "aws" "ali" "sugon"  |
| grouped | false | bool | weather group same task in to one job group. |