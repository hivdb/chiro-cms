# Chiro CMS

Content Management System of Chiro Website


## Folder Structure

- `pages/`: where we stores YAML files. A file starts with "page-" is
  automatically corresponded to a page on HIVDB website starts with "/page/".
  For example, "page-example.yml" is corresponded to page "/page/example/".
- `resources/`: for any YAML files in `pages/`, the user can specify a magic
  keyword `_resources` to refer to a file in this folder. The CMS system
  currently support both Markdown (`.md`) and HTML (`.html`) files. However
  please always choose Markdown as your first option since it's simpler and
  formatted more well than HTML file.
- `images/`: image files. You can use `$$CMS_PREFIX$$images/example.png` to
  refer to a image from this folder.
- `downloads/`: other files needed to be made downloadable. You can use
  `$$CMS_PREFIX$$downloads/example.txt` to refer to a file from this folder.

## Usage

### Prerequisites

- Docker: Follow https://docs.docker.com/install/ to install Docker on your
  computer.
- AWS IAM credential: Follow https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/setup-credentials.html
  to configure the credential. Your IAM must have the access to AWS S3.

If you have sudo access to the internal server `zhangfei`, it is already
configured under username "philip":

```bash
ssh zhangfei
sudo su - philip
cd ~/gitrepo/chiro-cms
```

### Deployment

Several deployment commands are available for different targets:

```bash
# build and deploy to development environment:
make deploy-dev

# build and deploy to covdb.stanford.edu:
make deploy-prod

# build and deploy to all targets:
make deploy-all
```

## Docker image

Chiro CMS uses the [same Docker image][docker-hub-link] used by HIVDB CMS.

[docker-hub-link]: https://hub.docker.com/r/hivdb/hivdb-cms-builder
